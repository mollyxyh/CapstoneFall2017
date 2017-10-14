import xlrd
import numpy as np

import sys
sys.path.append('../01 Pricing')
from SABR import SABR_model

def fitting(market_data, method='Hagan'):
    '''
    This function calibrates different versions of SABR models with the input data and does some preparation work for later over specification analysis, e.g. calibrating another three parameters with one fixed to specific values. The output results are saved in '05 Outputs' folder and '02 Fitter/parameters' folder.
    @var market_data: the input data to which SABR model needs to be fitted
    @val method: the specific version of SABR model we are using
    '''
    
    ######## inputs and outputs #########################################
    outvol = open('../05 Outputs/outvol_%s.csv' % (method), 'w')  # file output of volatilities
    vol_diff = open('../05 Outputs/vol_differences_%s.csv' % (method), 'w')  # file output differences between SABR and Market volatilities
    parameters = open('../05 Outputs/parameters_%s.csv' % (method), 'w')  # file output parameters

    file_input=xlrd.open_workbook('../04 Inputs/'+market_data)  # load market data
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

    outvol = open('../05 Outputs/outvol_%s.csv' % (method), 'a')  # file output of volatilities
    vol_diff = open('../05 Outputs/vol_differences_%s.csv' % (method),
                    'a')  # file output differences between SABR and Market volatilities
    parameters = open('../05 Outputs/parameters_%s.csv' % (method), 'a')  # file output parameters

    ######## Calibration ###################################################
    # Fitting SABR model with market_data
    sabr = SABR_model(label_ten, label_exp, num_strikes, label_strikes, strike_spreads, outvol, vol_diff, parameters)
    calibrates = sabr.calibration(starting_guess, F, K, expiries, MKT, 'auto', 'auto', method)
    print 'Calibrated results of input data:\n'
    sabr.SABR_vol_matrix(calibrates['alpha'], calibrates['beta'], calibrates['rho'], calibrates['nu'], F, K, expiries,
                         MKT, method)

    # Test Beta for over-specification analysis
    print '\nOver specification analysis with beta fixed:'
    for fix in [0, 0.3, 0.5, 0.7, 1]:
        calibrates = sabr.calibration(starting_guess, F, K, expiries, MKT, 1, fix, method)
        sabr.SABR_vol_matrix(calibrates['alpha'], calibrates['beta'], calibrates['rho'], calibrates['nu'], F, K,
                             expiries, MKT, method)

    # Test Rho for over-specification analysis
    print '\nOver specification analysis with rho fixed:'
    for fix in [0, -0.3, -0.5, -0.7, -0.9]:
        calibrates = sabr.calibration(starting_guess, F, K, expiries, MKT, 2, fix, method)
        sabr.SABR_vol_matrix(calibrates['alpha'], calibrates['beta'], calibrates['rho'], calibrates['nu'], F, K,
                             expiries, MKT, method)

    ## Test Nu
    # for fix in [0,0.3,0.5,0.7,1]:
    # calibrates=sabr.calibration(starting_guess,F,K,expiries,MKT,3,fix)
    # sabr.SABR_vol_matrix(calibrates['alpha'],calibrates['beta'],calibrates['rho'],calibrates['nu'],F,K,expiries,MKT)

    outvol.close()
    vol_diff.close()
    parameters.close()
