import xlrd
import numpy as np
import Pricing.SABR as SABR
from Pricing.Data_processor import data_reader,set_label,start_params

def fitting(market_data, file_name, method='Hagan',mark='auto'):
    '''
    This function calibrates different versions of SABR models with the input data and does some preparation work for later over specification analysis, e.g. calibrating another three parameters with one fixed to specific value. The output variable results is a dictionary that stores the ivols predicted, ivol gap from market ivols and parameters of the SABR model.
    @var market_data: the input data to which SABR model needs to be fitted
    @var method: the specific version of SABR model we are using
    @var mark: representing if equal constraint is applied and of which parameter the constraint is applied
    '''
    
    ######## inputs and outputs #########################################
    data = data_reader(market_data,file_name)

    # set starting parameters
    starting_guess = np.array([0.001, 0.5, 0, 0.001])
    start_params(data['F'],starting_guess)

    ######## set labels ###################################################
    labels = set_label(data['strike_spreads'],data['num_strikes'],data['expiries'],data['tenors'])
    
    ######## Calibration ###################################################
    results={} # set output dictionary
    
    if mark=='auto':
        # Fit SABR model with market_data
        sabr=SABR.SABR_model(labels['label_ten'],labels['label_exp'],data['num_strikes'],labels['label_strikes'],data['strike_spreads'])
        calibrates_auto=sabr.calibration(starting_guess, data['F'], data['K'], data['expiries'], data['MKT'], 'auto', 'auto', method)
        results['auto']=sabr.SABR_vol_matrix(calibrates_auto['alpha'], calibrates_auto['beta'], calibrates_auto['rho'], calibrates_auto['nu'], data['F'], data['K'], data['expiries'], data['MKT'], method)
        results['jacmat']=calibrates_auto['jacmat'] # stored for collinearity analysis
        
    elif mark=='beta':
        # Test Beta for over-specification analysis
        for fix in [0, 0.3, 0.5, 0.7, 1]:
            sabr=SABR.SABR_model(labels['label_ten'],labels['label_exp'],data['num_strikes'],
                                 labels['label_strikes'],data['strike_spreads'])
            calibrates=sabr.calibration(starting_guess, data['F'], data['K'], data['expiries'], data['MKT'], 1, fix, method)
            results[fix]=sabr.SABR_vol_matrix(calibrates['alpha'], calibrates['beta'], calibrates['rho'], calibrates['nu'], data['F'], data['K'], data['expiries'], data['MKT'], method)
    
    elif mark=='rho':
        # Test Rho for over-specification analysis
        for fix in [0, -0.3, -0.5, -0.7, -0.9]:
            sabr=SABR.SABR_model(labels['label_ten'],labels['label_exp'],data['num_strikes'],
                                 labels['label_strikes'],data['strike_spreads'])
            calibrates = sabr.calibration(starting_guess, data['F'], data['K'], data['expiries'], data['MKT'], 2, fix, method)
            results[fix]=sabr.SABR_vol_matrix(calibrates['alpha'], calibrates['beta'], calibrates['rho'], calibrates['nu'], data['F'], data['K'], data['expiries'], data['MKT'], method)
    
    elif mark=='alpha':   
        # Test Alpha for over-specification analysis
        for fix in [0.2, 0.4, 0.6]:
            sabr=SABR.SABR_model(labels['label_ten'],labels['label_exp'],data['num_strikes'],
                                 labels['label_strikes'],data['strike_spreads'])
            calibrates = sabr.calibration(starting_guess, data['F'], data['K'], data['expiries'], data['MKT'], 0, fix, method)
            results[fix]=sabr.SABR_vol_matrix(calibrates['alpha'], calibrates['beta'], calibrates['rho'], calibrates['nu'], data['F'], data['K'], data['expiries'], data['MKT'], method)
    
    elif mark=='vega':
        # Test vega for over-specification analysis
        for fix in [0.2,0.4,0.6]:
            sabr=SABR.SABR_model(labels['label_ten'],labels['label_exp'],data['num_strikes'],
                                 labels['label_strikes'],data['strike_spreads'])
            calibrates=sabr.calibration(starting_guess,data['F'], data['K'], data['expiries'], data['MKT'],3, fix, method)
            results[fix]=sabr.SABR_vol_matrix(calibrates['alpha'], calibrates['beta'], calibrates['rho'], calibrates['nu'], data['F'], data['K'], data['expiries'], data['MKT'], method)
    
    return results

