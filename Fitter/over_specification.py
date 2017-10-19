import numpy as np
from matplotlib import pyplot as plt

def over_specification(check):
    '''
    This function plots ivols smile by Hagan Lognormal SABR with different equal constraints and the market ivols smile.
    @var check: a dictionary returned by fitting.py, which stores the calibration results with a specific parameter fixed to specific values.
    '''
    keys=check.keys() #a list of values to which the parameter is fixed
    K_spreads=check[keys[0]]['ivols'].columns.values.tolist()[3:] #set K spreads as xlabels
    K_spreads[K_spreads.index('ATM')]=0
    for key in keys:
        results=check[float(key)] #calibration result with specific equal constraint
        vols=results['ivols'].loc[27,:].values[3:] #ivols predicted for different K spreads
        vols_diff=results['ivols_diff'].loc[27,:].values[3:] #ivols gaps for different K spreads
        MKT=vols-vols_diff #market ivols for different K spreads
        plt.plot(K_spreads, vols*100, linestyle='solid',label='parameter='+str(key)) #calibrated implied vols by SABR
    plt.plot(K_spreads, MKT*100, linestyle='--', color='black', label='Market') #market implied vols
    plt.title('Over-specification Analysis')
    plt.xlabel('strike spreads (bps)')
    plt.ylabel('implied vols (%)')
    plt.legend()
    plt.show()
    