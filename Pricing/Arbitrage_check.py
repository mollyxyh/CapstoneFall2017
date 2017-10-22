import numpy as np
import pandas as pd
import xlrd
from matplotlib import pyplot as plt
import Pricing.black_pricing as BS
import Pricing.SABR as SABR
from Pricing.Data_processor import data_reader,set_label,start_params
from Pricing.pdf_calculation import numerical_pdf
from Fitter.fitting import fitting,objfunc,calibration

def arbitrage_check(market_data,file_name):
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

