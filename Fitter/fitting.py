import xlrd
import numpy as np
import Pricing.SABR as SABR

def fitting(market_data, method='Hagan'):
    '''
    This function calibrates different versions of SABR models with the input data and does some preparation work for later over specification analysis, e.g. calibrating another three parameters with one fixed to specific values. The output results are saved in '05 Outputs' folder and '02 Fitter/parameters' folder.
    @var market_data: the input data to which SABR model needs to be fitted
    @val method: the specific version of SABR model we are using
    '''
    
    ######## inputs and outputs #########################################
    file_input=xlrd.open_workbook('../Inputs/'+market_data)  # load market data
    Market_data=file_input.sheet_by_name('Swaptions data')  # file input forward rates

    ######## set swaptions characteristics ###############################
    strike_spreads=[]
    j = 0
    while True:
        try:
            strike_spreads.append(int(Market_data.cell(1, 3 + j).value))
            j = j + 1
        except:
            break
    num_strikes = len(strike_spreads)

    expiries = []
    i = 0
    while True:
        try:
            expiries.append(Market_data.cell(2 + i, 1).value)
            i = i + 1
        except:
            break

    tenors = []
    i = 0
    while True:
        try:
            tenors.append(Market_data.cell(2 + i, 0).value)
            i = i + 1
        except:
            break

    # to create the ATM forward rates
    F = []
    i = 0
    while True:
        try:
            F.append(Market_data.cell(2 + i, 2).value)
            i = i + 1
        except:
            break

    # to create the strike grid
    K = np.zeros((len(F), num_strikes))
    for i in range(len(F)):
        for j in range(num_strikes):
            K[i][j] = F[i] + 0.0001 * (strike_spreads[j])

            # to create market volatilities
    MKT = np.zeros((len(F), num_strikes))
    for i in range(len(F)):
        for j in range(num_strikes):
            MKT[i][j] = Market_data.cell(2 + i, 3 + j).value

    # set starting parameters
    global alpha, beta, rho, nu, jacmat
    starting_guess = np.array([0.001, 0.5, 0, 0.001])
    alpha = len(F) * [starting_guess[0]]
    beta = len(F) * [starting_guess[1]]
    rho = len(F) * [starting_guess[2]]
    nu = len(F) * [starting_guess[3]]
    jacmat = len(F) * [starting_guess[3]]

    ######## set labels ###################################################
    exp_dates = len(expiries) * [0]
    for i in range(len(expiries)):
        if expiries[i] < 1:
            exp_dates[i] = str(int(round(12 * expiries[i]))) + 'm'
        else:
            exp_dates[i] = str(int(round(expiries[i]))) + 'y'
            if expiries[i] - round(expiries[i]) > 0:
                exp_dates[i] = exp_dates[i] + str(int(round((12 * (round(expiries[i], 2) - int(expiries[i])))))) + 'm'
            elif expiries[i] - round(expiries[i]) < 0:
                exp_dates[i] = str(int(round(tenors[i])) - 1) + 'y'
                exp_dates[i] = exp_dates[i] + str(int(round((12 * (round(expiries[i], 2) - int(expiries[i])))))) + 'm'

    ten_dates = len(tenors) * [0]
    for i in range(len(tenors)):
        if tenors[i] < 1:
            ten_dates[i] = str(int(round(12 * tenors[i]))) + 'm'
        else:
            ten_dates[i] = str(int(round(tenors[i]))) + 'y'
            if tenors[i] - round(tenors[i]) > 0:
                ten_dates[i] = ten_dates[i] + str(int(round((12 * (round(tenors[i], 2) - int(tenors[i])))))) + 'm'
            elif tenors[i] - round(tenors[i]) < 0:
                ten_dates[i] = str(int(round(tenors[i])) - 1) + 'y'
                ten_dates[i] = ten_dates[i] + str(int(round((12 * (round(tenors[i], 2) - int(tenors[i])))))) + 'm'

    label_exp = exp_dates
    label_ten = ten_dates
    label_strikes = num_strikes * [0]
    for i in range(num_strikes):
        if strike_spreads[i] == 0:
            label_strikes[i] = 'ATM'
        else:
            label_strikes[i] = str(strike_spreads[i])

    ######## Calibration ###################################################
    # Fitting SABR model with market_data
    sabr=SABR.SABR_model(label_ten,label_exp,num_strikes,label_strikes,strike_spreads)
    calibrates_auto=sabr.calibration(starting_guess, F, K, expiries, MKT, 'auto', 'auto', method)
    auto=sabr.SABR_vol_matrix(calibrates_auto['alpha'], calibrates_auto['beta'], calibrates_auto['rho'], calibrates_auto['nu'], F, K, expiries,
                         MKT, method)

    # Test Beta for over-specification analysis
    beta_check={}
    for fix in [0, 0.3, 0.5, 0.7, 1]:
        sabr=SABR.SABR_model(label_ten,label_exp,num_strikes,label_strikes,strike_spreads)
        calibrates_beta=sabr.calibration(starting_guess, F, K, expiries, MKT, 1, fix, method)
        beta_check[fix]=sabr.SABR_vol_matrix(calibrates_beta['alpha'], calibrates_beta['beta'], calibrates_beta['rho'], calibrates_beta['nu'], F, K,
                             expiries, MKT, method)

    # Test Rho for over-specification analysis
    rho_check={}
    for fix in [0, -0.3, -0.5, -0.7, -0.9]:
        sabr=SABR.SABR_model(label_ten,label_exp,num_strikes,label_strikes,strike_spreads)
        calibrates = sabr.calibration(starting_guess, F, K, expiries, MKT, 2, fix, method)
        rho_check[fix]=sabr.SABR_vol_matrix(calibrates['alpha'], calibrates['beta'], calibrates['rho'], calibrates['nu'], F, K,
                             expiries, MKT, method)
        
    # Test Alpha for over-specification analysis
    alpha_check={}
    for fix in [0.2, 0.4, 0.6]:
        sabr=SABR.SABR_model(label_ten,label_exp,num_strikes,label_strikes,strike_spreads)
        calibrates = sabr.calibration(starting_guess, F, K, expiries, MKT, 0, fix, method)
        alpha_check[fix]=sabr.SABR_vol_matrix(calibrates['alpha'], calibrates['beta'], calibrates['rho'], calibrates['nu'], F, K,
                             expiries, MKT, method)

    # Test Nu for over-specification analysis
    vega_check={}
    for fix in [0.2,0.4,0.6]:
        sabr=SABR.SABR_model(label_ten,label_exp,num_strikes,label_strikes,strike_spreads)
        calibrates=sabr.calibration(starting_guess,F,K,expiries,MKT,3, fix, method)
        vega_check[fix]=sabr.SABR_vol_matrix(calibrates['alpha'], calibrates['beta'], calibrates['rho'], calibrates['nu'], F, K,
                             expiries, MKT, method)
    return {'auto':auto,'beta_check':beta_check,'rho_check':rho_check,'alpha_check':alpha_check,'vega_check':vega_check,'jacmat':calibrates_auto['jacmat']}
        