import numpy as np
import yaml

DEFAULT_ORDER = "AHCOMLNP"

def make_monod_function(c_half) -> callable:
    if c_half is None:
        return lambda c : 1.0 if c > 0 else 0.0
    else: 
        return lambda c : max( 0, ( c / ( c_half + c ) ) )

def gen_ivp_terminate_event(idx):
    def event(t,y):
        return y[idx]
    
    event.terminal = True
    event.direction = -1
    return event

IVP_SOLUTION_TERMINATION_EVENTS = [ gen_ivp_terminate_event(i) for i,n in enumerate(DEFAULT_ORDER) if n in ["N","P"] ] # termination if N or P run out

class suspended_cultures_ode:
    def __init__(self):
        ### initialize dictionaries for parameter groups
        self.monods = {}
        self.base_growth_rates = {}
        self.inverse_yields = {}
        self.alpha = {}
        self.beta = {}
        self.exchange_coeff = {}
        self.equilibrium_gas_concentration = {}

    def load_reaction_params(self, constants_filename:str, autotrophs:str="A_platensis", heterotrophs:str="M_capsulatus") -> bool:
        CONSTANTS = None
        with open(constants_filename, 'r') as handle: CONSTANTS = yaml.safe_load(handle)
        if CONSTANTS is None: return False

        for b1,b2 in dict(A=autotrophs, H=heterotrophs).items():
            ### maximum growth rates
            self.base_growth_rates[b1] = np.mean( np.log(2) / np.array(
                    [ arr[0] for arr in CONSTANTS["doubling_times"][b2] ]
            ) )

            for n in ["C","O","M","L","N","P"]:
                ### Monod function components
                if n in CONSTANTS["half_saturation_concentrations"][b2]:
                    c_half = CONSTANTS["half_saturation_concentrations"][b2][n]
                    if type(c_half) is not float:
                        print(f"warning: type of c_half is {type(c_half)}, not float.")
                        c_half = float(c_half)
                    self.monods[f"{b1}<-{n}"] = make_monod_function(c_half)
                
                ### inverse yields
                if n in CONSTANTS["inverse_yields"][b2]:
                    self.inverse_yields[f"{b1}<-{n}"] = CONSTANTS["inverse_yields"][b2][n]

                ### gas production coefficients
                if n in CONSTANTS["gas_prod_coeff_growth_assoc"][b2]:
                    self.alpha[f"{b1}->{n}"] = CONSTANTS["gas_prod_coeff_growth_assoc"][b2][n]
                if n in CONSTANTS["gas_prod_coeff_non_growth_assoc"][b2]:
                    self.beta[f"{b1}->{n}"] = CONSTANTS["gas_prod_coeff_non_growth_assoc"][b2][n]
        del b1,b2,n

        return True
    
    def load_gas_exchange_params(self, constants_filename:str, params_filename:str) -> bool:
        CONSTANTS = None
        with open(constants_filename, 'r') as handle: CONSTANTS = yaml.safe_load(handle)
        if CONSTANTS is None: return False

        PARAMS = None
        with open(params_filename, 'r') as handle: PARAMS = yaml.safe_load(handle)
        if PARAMS is None: return False

        volume = 1e-3*PARAMS["total_volume_litres"]
        interface = 1e-4*PARAMS["surface_area_sqcm"]
        supply_composition = { n:v/(1+PARAMS["methane_ratio_in_supply"]) for n,v in CONSTANTS["standard_atmosphere_composition"].items() }
        supply_composition["M"] = PARAMS["methane_ratio_in_supply"]
        
        for n in supply_composition.keys():
            partial_pressure = 100*supply_composition[n]*PARAMS["gas_supply_pressure_hpa"]
            henry_law_coeff = np.mean( [ arr[0] for arr in CONSTANTS["henry_law_coeffs"][n] ] )
            self.equilibrium_gas_concentration[n] = henry_law_coeff * partial_pressure

            diffusion_coeff = np.mean( [ arr[0] for arr in CONSTANTS["diffusion_coeffs_in_water"][n] ] )
            stirring_time = PARAMS["stirring_time_scale_seconds"]
            mass_transf_coeff = 2*np.sqrt( diffusion_coeff / np.pi / stirring_time )

            self.exchange_coeff[n] = mass_transf_coeff * interface / volume
        del n

        return True
        
    def _mu(self, A,H,C,O,M,L,N,P):
        mu_A = self.base_growth_rates["A"]
        mu_H = self.base_growth_rates["H"]
        
        ### multiplicative model
        for key,value in dict(C=C, L=L, N=N, P=P).items():
            mu_A *= self.monods[f"A<-{key}"]( value )
        for key,value in dict(O=O, M=M, N=N, P=P).items():
            mu_H *= self.monods[f"H<-{key}"]( value )
        
        return mu_A, mu_H
    
    def _consumption_rates(self, dA, dH):
        ### initialize components of the output vector
        dC,dO,dM,dL,dN,dP = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0

        ### uptake by autotrophs
        dC -= dA * self.inverse_yields["A<-C"] ### CO2
        ### TODO: light
        dN -= dA * self.inverse_yields["A<-N"] ### nitrogen
        dP -= dA * self.inverse_yields["A<-P"] ### phosphorous
        
        ### uptake by heterotrophs
        dO -= dH * self.inverse_yields["H<-O"] ### oxygen
        dM -= dH * self.inverse_yields["H<-M"] ### methane
        dN -= dH * self.inverse_yields["H<-N"] ### nitrogen
        dP -= dH * self.inverse_yields["H<-P"] ### phosphorous
        
        return dC,dO,dM,dL,dN,dP
    
    def _gas_production(self, A, H, dA, dH):
        dO = dA * self.alpha["A->O"] + A * self.beta["A->O"]
        dC = dH * self.alpha["H->C"] + H * self.beta["H->C"]
        return dC, dO
    
    def _gas_exchange(self, C, O, M, N):
        dC = self.exchange_coeff["C"] * ( self.equilibrium_gas_concentration["C"] - C )
        dO = self.exchange_coeff["O"] * ( self.equilibrium_gas_concentration["O"] - O )
        dM = self.exchange_coeff["M"] * ( self.equilibrium_gas_concentration["M"] - M )
        #dN = self.exchange_coeff["N"] * ( self.equilibrium_gas_concentration["N"] - N )
        dN = 0.0
        #ignore nitrogen because it's not supplied in elemental gas form
        return dC, dO, dM, dN

    ### to calculate the derivatives
    def __call__(self, t, state):
        ### unpack state variables
        A,H,C,O,M,L,N,P = state
        
        ### calculate growth rates
        mu_A, mu_H = self._mu(*state)

        ### proliferation
        dA = mu_A * A
        dH = mu_H * H

        ### nutrient removal
        dC,dO,dM,dL,dN,dP = self._consumption_rates(dA, dH)

        ### gas production
        prod_dC, prod_dO = self._gas_production(A,H,dA,dH)
        dC += prod_dC
        dO += prod_dO

        ### gas exchange
        exc_dC, exc_dO, exc_dM, exc_dN = self._gas_exchange(C, O, M, N)
        dC += exc_dC
        dO += exc_dO
        dM += exc_dM
        dN += exc_dN

        ### return derivative vector
        return [ dA, dH, dC, dO, dM, dL, dN, dP ]

def make_initial_state(ymlfile:str):
    with open(ymlfile, 'r') as ymlfile:
        PARAMS = yaml.safe_load(ymlfile)
        return [ PARAMS[key] for key in "AHCOMLNP" ]

def update_state_with_gas_at_equilibrium(state, ode_system:suspended_cultures_ode):
    res = [ s for s in state ]
    for gas in "COM": # ignore nitrogen
        idx = { gas:i for i,gas in enumerate(DEFAULT_ORDER) }[gas]
        res[idx] = ode_system.equilibrium_gas_concentration[gas]
    return res


### Tests
def test():
    pass

if __name__ == "__main__":
    test()
