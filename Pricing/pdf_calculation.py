from scipy.stats import lognorm
import numpy as np
import pandas as pd
import xlrd
import math
from matplotlib import pyplot as plt
import Pricing.black_pricing as BS
import Pricing.SABR as SABR
from Pricing.Data_processor import data_reader,set_label,start_params
from Fitter.fitting import fitting,objfunc,calibration

def numerical_pdf(market_data,file_name,h=0.0001,vol_method='Hagan'):
    data = data_reader(market_data,file_name)    #read data
    
    # generate new strike spreads and strike prices with smaller intervals for approximation
    new_strike_spreads = np.arange(data['strike_spreads'][0],data['strike_spreads'][-1]+1,1)
    new_num_strikes = len(new_strike_spreads)
    new_labels = set_label(new_strike_spreads,new_num_strikes,data['expiries'],data['tenors'])  #set new labels
    
    new_K = np.zeros((len(data['F']), new_num_strikes))  #generate new strike prices
    for i in range(len(data['F'])):
        for j in range(new_num_strikes):
            new_K[i][j] = data['F'][i] + 0.0001 * (new_strike_spreads[j])
            
    label_strikes=['tenor','expiry','F']
    for j in range(len(new_labels['label_strikes'])):
        label_strikes.append(new_labels['label_strikes'][j])
    
    # caliberate alpha,beta,pho,nu with market_data for computing SABR vols
    starting_guess = np.array([0.001, 0.5, 0, 0.001])
    start_params(data['F'],starting_guess)
    
    calibrates_auto=calibration(starting_guess, data['F'], data['K'], data['expiries'], data['MKT'], 'auto', 'auto', vol_method)
    
    # computing SABR vols for new strike prices
    sabr=SABR.SABR_model(new_labels['label_ten'],new_labels['label_exp'],new_num_strikes,
                         new_labels['label_strikes'],new_strike_spreads)
    new_ivols = sabr.SABR_vol_matrix(calibrates_auto['alpha'], calibrates_auto['beta'], calibrates_auto['rho'], 
                            calibrates_auto['nu'], data['F'], new_K, data['expiries'], vol_method)

    # pricing for different options
    bs = BS.BS_pricing(new_labels['label_ten'],new_labels['label_exp'],new_num_strikes,new_labels['label_strikes'],new_strike_spreads)
    price = bs.bs_matrix(data['F'], new_K, data['expiries'],new_ivols,1,0)
    
    # approximate second derivative with 2nd order central method
    numerical_pdf = []
    for i in range(len(data['F'])):
        #p = price.loc[i][3:]
        temp=[new_labels['label_ten'][i],new_labels['label_exp'][i],data['F'][i]]
        for j in range(3,len(label_strikes)):
            if j==3 or j==len(label_strikes)-1:
                pdf = 0
            else:
                pdf = (price.loc[i].values[j+1]-2*price.loc[i].values[j]+price.loc[i].values[j-1])/h**2
            temp.append(pdf)
        numerical_pdf.append(temp)
    pdf_ = pd.DataFrame(data=numerical_pdf,columns=label_strikes)
    result = {'pdf':pdf_,'strike_price':new_K}
    
    K = result['strike_price']
    pdf = result['pdf']
    tenor = pdf['tenor']
    expiry = pdf['expiry']

    fig, ax = plt.subplots()
    for i in range(pdf.shape[0]):
        ax.plot(K[i],pdf.loc[i][3:],label=(tenor[i]+'; '+expiry[i]))
        ax.legend()
    plt.xlabel('strike price')
    plt.ylabel('lognormal pdf')
    plt.title('Arbitrage Check')
    plt.show()
    return result

def lognormal_pdf(x,sigma,drift):
    result = 1/(x*sigma*math.sqrt(2*pi)) * exp*(-(math.log(x)-drift)**2/(2*sigma*sigma))
    return result