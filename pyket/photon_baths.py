
def lorentzian(x,x0,width):
    x = (x-x0)/(width/2)
    return _standardized_lorentzian(x)
    

def _standardized_lorentzian(x):

    return 1/(x+x**2)
