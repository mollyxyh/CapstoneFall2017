import numpy as np
from matplotlib import pyplot as plt
from Fitter.fitter import Fitter
from Pricing.SABR import SABR_model

def over_specification(fitter,check):
    '''
    This function plots ivols smile by Hagan Lognormal SABR with different equal constraints and the market ivols smile.
    @var check: a dictionary returned by fitter.py, which stores the calibration results with a specific parameter fixed to specific values.
    '''
    keys=check.keys() #a list of values to which the parameter is fixed
    #K_spreads=[-150,-100,-50,-25,0,25,50,100,150]
    for key in keys:
        results=check[key] #calibration result with specific equal constraint
        ivols=results.iloc[10,3:].values.tolist()
        plt.plot(ivols,linestyle='solid',label='parameter='+str(key)) #calibrated implied vols by SABR
    plt.title('Over-specification Analysis')
    plt.xlabel('K spreads (bps)')
    plt.ylabel('ivols (%)')
    plt.legend()
    plt.show()
    
def get_check(fitter,fix_par,fix_no):
    check={}
    for no in fix_no:
        results=fitter.calibration(eqc=[fix_par,no])
        check[no]=fitter.ivol_SABR(results['alpha'],results['beta'],results['rho'],results['nu'])
    return check