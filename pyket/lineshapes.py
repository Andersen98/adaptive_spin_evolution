
def lorentzian(x,x0,width):
    assert(width > 0)
    x = (x-x0)/(width/2)
    return _standardized_lorentzian(x)
    

def _standardized_lorentzian(x):

    return 1/(1+x**2)
