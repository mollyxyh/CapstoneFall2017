import numpy as np
from matplotlib import pyplot as plt
from Fitter.fitter import Fitter
from Pricing.SABR import SABR_model

def over_specification(fitter,check):
    '''
    This function plots ivols smile by Hagan Lognormal SABR with different equal constraints and the market ivols smile.
    @var fitter: a fitter object
    @var check: a dictionary returned by fitter.py, which stores the calibration results with a specific parameter fixed to specific values.
    '''
    keys=check.keys()
    K_spread=fitter.K_spread #x label
    MKT=fitter.MKT
    ivols_mkt=MKT[10]
    ivols_mkt=[item*100 for item in ivols_mkt]
    for key in keys:
        results=check[key] #calibration result with specific equal constraint
        ivols=results.iloc[10,3:].values.tolist() #take 5Y tenor and 1Y expiry as an example
        ivols=[item*100 for item in ivols] #y label
        plt.plot(K_spread,ivols,linestyle='solid',label='parameter='+str(key))
    plt.plot(K_spread,ivols_mkt,linestyle='--',color='black',label='MKT') #market ivol
    plt.title('Over-specification Analysis')
    plt.xlabel('K spreads (bps)')
    plt.ylabel('ivols (%)')
    plt.legend()
    plt.show()
    
def get_check(fitter,fix_par,fix_no,method='Hagan'):
    '''
    This function generates a dictionary composed of ivols dataframes predicted under SABR model with different calibrations.
    @var fitter: a fitter object
    @var fix_par: parameter that needs to be fixed when calibrated
    @var fix_no: a vector of numbers to which fix_par needs to be fixed
    '''
    check={}
    for no in fix_no:
        results=fitter.calibration(eqc=[fix_par,no],method=method)
        check[no]=fitter.ivol_SABR(results['alpha'],results['beta'],results['rho'],results['nu'],method=method)
    return check