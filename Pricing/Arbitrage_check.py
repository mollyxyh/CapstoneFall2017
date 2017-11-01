import numpy as np
import pandas as pd
import xlrd
from matplotlib import pyplot as plt
#from Pricing.black_pricing import BSPricer_SABR
#from Pricing.SABR import SABR_model
#from Pricing.Data_processor import data_reader,set_label,start_params
from Pricing.pdf_calculation import numerical_pdf

def arbitrage(alpha,F,K,expiry,r=0,isCall=1,h=0.0001,vol_method='Hagan',vol_dist='lognormal'):
    ab = []
    temp = []
    pdf = numerical_pdf(alpha,F,K,expiry,isCall,r,h,vol_method,vol_dist)
    print('SABR volatility method:',vol_method)
    print('volatility distribution:',vol_dist)
    for i in range(pdf.shape[0]):
        for j in range(pdf.shape[1]):
            if pdf[i][j] <= 0:
                x = 1
                print('Expiry=',expiry[i],'Strike price=',K[i][j],'p.d.f=',pdf[i][j],'Yes')
            else:
                x = 0
                print('Expiry=',expiry[i],'Strike price=',K[i][j],'p.d.f=',pdf[i][j],'No')
            temp.append(x)
        ab.append(temp)
    result = np.array(ab)
    return result

"""
    def pdf_plot():
        fig, ax = plt.subplots()
        for i in range(pdf.shape[0]):
            ax.plot(K[i],pdf.loc[i][3:],label=(tenor[i]+'; '+expiry[i]))
            ax.legend()
        plt.xlabel('strike price')
        plt.ylabel('lognormal pdf')
        plt.title('Arbitrage Check')
        plt.show()
        return result
"""
    

        
    
"""
def arbitrage_check(K,pdf):
    data = data_reader(market_data,file_name)    #read data
    
    new_strike_spreads = np.arange(data['strike_spreads'][0],data['strike_spreads'][-1]+1,1)
    new_num_strikes = len(new_strike_spreads)
    new_labels = set_label(new_strike_spreads,new_num_strikes,data['expiries'],data['tenors'])  #set new labels
    
    label_strikes=['tenor','expiry','F']
    for j in range(len(new_labels['label_strikes'])):
        label_strikes.append(new_labels['label_strikes'][j])
        
    arbitrage = []
    results = numerical_pdf(market_data,file_name,h=0.0001,vol_method='Hagan')
    pdf = results['pdf']
    for i in range(pdf.shape[0]):
        temp=[new_labels['label_ten'][i],new_labels['label_exp'][i],data['F'][i]]
        for j in range(3,pdf.shape[1]):
            if j==3 or j==pdf.shape[1]-1:
                x = 'No data'
            else:
                if pdf.loc[i].values[j] < 0 or pdf.loc[i].values[j]==0:
                    x = 1  # Arbitrage violation
                else:
                    x = 0  # No arbitrage
            temp.append(x)
        arbitrage.append(temp)
    results = pd.DataFrame(data=arbitrage,columns=label_strikes)
    return results
"""

