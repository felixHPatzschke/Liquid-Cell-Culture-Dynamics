import numpy as np

# needed to package into proper function for some weird late binding reason
def make_monod_function(c_half) -> callable:
    if c_half is None:
        return lambda c : 1.0 if c > 0 else 0.0
    else: 
        return lambda c : max( 0, ( c / ( c_half + c ) ) )
    
def make_multiplicative_monod_function(  ):
    pass
