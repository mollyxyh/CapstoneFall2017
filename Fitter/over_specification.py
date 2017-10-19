import numpy as np
from matplotlib import pyplot as plt

def over_specification(check):
    '''
    This function plots implied vols under Hagan Lognormal SABR with parameter[alpha,beta,rho,nu] fixed to specific values vs. market implied vols in one plot.
    @var check: result dictionary outputed by fitting.py, which stores the calibration results with specific parameter fixed to specific values.
    '''
    keys=check.keys()
    K_spreads=check[keys[0]]['ivols'].columns.values.tolist()[3:] #set K spreads as xlabels
    K_spreads[K_spreads.index('ATM')]=0
    for key in keys:
        results=check[float(key)]
        vols=results['ivols'].loc[27,:].values[3:]
        vols_diff=results['ivols_diff'].loc[27,:].values[3:]
        MKT=vols-vols_diff
        plt.plot(K_spreads, vols*100, linestyle='solid',label='parameter='+str(key)) #calibrated implied vols
    plt.plot(K_spreads, MKT * 100, linestyle='--', color='black', label='Market') #market implied vols
    plt.title('Over-specification Analysis')
    plt.xlabel('strike spreads (bps)')
    plt.ylabel('implied vols (%)')
    plt.legend()
    plt.show()
    