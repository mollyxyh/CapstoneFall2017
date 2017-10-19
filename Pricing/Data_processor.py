import xlrd

def data_reader(file_name,sheet_name):
    file_input=xlrd.open_workbook('../Inputs/'+file_name)  # import market data
    Market_data=file_input.sheet_by_name(sheet_name)  # file input forward rates
    
     ######## set characteristics ###############################
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


def set_label():
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
