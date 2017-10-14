import numpy as np
from matplotlib import pyplot as plt

def check_beta(vols,vol_diff):
    '''
    This function plots implied vols under Hagan Lognormal SABR with beta fixed to [0,0.3,0.5,0.7,1] vs. market implied vols.
    @var vols: calibrated vol data from SABR_calibration.ipynb
    @val vol_diff: calibrated vol difference data (vol-mkt) from SABR_calibration.ipynb
    '''
    K=[-150,-100,-50,-25,0,25,50,100,150]
    bs = ['auto', 0, 0.3, 0.5, 0.7, 1] # auto means optimization without equality constraints, which means optimal calibration
    for b in [0, 1, 2, 3, 4, 5]:
        vols_tar = vols.loc[28 + 38 * b].values[0]
        vols_tar = vols_tar.split(';')[2:-1]
        vols_tar = np.array([float(item) for item in vols_tar])
        diff = vol_diff.loc[28 + 38 * b].values[0]
        diff = diff.split(';')[2:-1]
        diff = np.array([float(item) for item in diff])
        MKT = vols_tar - diff

        plt.plot(K, vols_tar * 100, linestyle='solid', marker='x', color='black', label='Hagan Lognormal SABR') #calibrated implied vols
        plt.plot(K, MKT * 100, linestyle='--', color='blue', label='Market') #market implied vols
        plt.title('Beta checking with beta=' + str(bs[b]))
        plt.xlabel('strike')
        plt.ylabel('implied vols (%)')
        plt.legend()
        plt.show()
        
def check_rho(vols,vol_diff):
    '''
    This function plots implied vols under Hagan Lognormal SABR with rho fixed to [0,0.3,0.5,0.7,0.9] vs. market implied vols.
    @var vols: calibrated vol data from SABR_calibration.ipynb
    @val vol_diff: calibrated vol difference data (vol-mkt) from SABR_calibration.ipynb
    '''
    K=[-150,-100,-50,-25,0,25,50,100,150]
    bs = [0, -0.3, -0.5, -0.7, -0.9]
    for b in [0, 1, 2, 3, 4]:
        vols_tar = vols.loc[28 + 38 * (b+6)].values[0]
        vols_tar = vols_tar.split(';')[2:-1]
        vols_tar = np.array([float(item) for item in vols_tar])
        diff = vol_diff.loc[28 + 38 * b].values[0]
        diff = diff.split(';')[2:-1]
        diff = np.array([float(item) for item in diff])
        MKT = vols_tar - diff

        plt.plot(K, vols_tar * 100, linestyle='solid', marker='x', color='black', label='Hagan Lognormal SABR')
        plt.plot(K, MKT * 100, linestyle='--', color='blue', label='Market')
        plt.title('Rho checking with rho=' + str(bs[b]))
        plt.xlabel('strike')
        plt.ylabel('implied vols (%)')
        plt.legend()
        plt.show()