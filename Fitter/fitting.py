import xlrd
import numpy as np
import pandas as pd
import math
import Pricing.SABR as SABR
from scipy.optimize import minimize
from Pricing.Data_processor import data_reader,set_label,start_params

def fitting(market_data, file_name, method='Hagan',mark='auto'):
    '''
    This function calibrates different versions of SABR models with the input data and does some preparation work for later
    over specification analysis, e.g. calibrating another three parameters with one fixed to specific value. The output
    variable results is a dictionary that stores the ivols predicted, ivol gap from market ivols and parameters of the SABR
    model.
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
        sabr=SABR.SABR_model(labels['label_ten'],labels['label_exp'],data['num_strikes'],
                             labels['label_strikes'],data['strike_spreads'])
        calibrates_auto=calibration(starting_guess, data['F'], data['K'], data['expiries'], 
                                    data['MKT'], 'auto', 'auto', method)
        results['auto_vols']=sabr.SABR_vol_matrix(calibrates_auto['alpha'], calibrates_auto['beta'], calibrates_auto['rho'],
                                                  calibrates_auto['nu'], data['F'], data['K'], data['expiries'], method)
        results['auto_params']=calibrates_auto['params']
        results['jacmat']=calibrates_auto['jacmat'] # stored for collinearity analysis
        
    elif mark=='beta':
        # Test Beta for over-specification analysis
        for fix in [0, 0.3, 0.5, 0.7, 1]:
            sabr=SABR.SABR_model(labels['label_ten'],labels['label_exp'],data['num_strikes'],
                                 labels['label_strikes'],data['strike_spreads'])
            calibrates=calibration(starting_guess, data['F'], data['K'], data['expiries'], data['MKT'], 1, fix, method)
            results[fix+'_vols']=sabr.SABR_vol_matrix(calibrates['alpha'], calibrates['beta'], calibrates['rho'], 
                                                      calibrates['nu'], data['F'], data['K'], data['expiries'], method)
            results[fix+'_params']=calibrates['params']
    
    elif mark=='rho':
        # Test Rho for over-specification analysis
        for fix in [0, -0.3, -0.5, -0.7, -0.9]:
            sabr=SABR.SABR_model(labels['label_ten'],labels['label_exp'],data['num_strikes'],
                                 labels['label_strikes'],data['strike_spreads'])
            calibrates = calibration(starting_guess, data['F'], data['K'], data['expiries'], data['MKT'], 2, fix, method)
            results[fix+'_vols']=sabr.SABR_vol_matrix(calibrates['alpha'], calibrates['beta'],calibrates['rho'],
                                                      calibrates['nu'],data['F'], data['K'], data['expiries'], method)
            results[fix+'_params']=calibrates['params']
    
    elif mark=='alpha':   
        # Test Alpha for over-specification analysis
        for fix in [0.2, 0.4, 0.6]:
            sabr=SABR.SABR_model(labels['label_ten'],labels['label_exp'],data['num_strikes'],
                                 labels['label_strikes'],data['strike_spreads'])
            calibrates = calibration(starting_guess, data['F'], data['K'], data['expiries'],data['MKT'], 0, fix, method)
            results[fix+'_vols']=sabr.SABR_vol_matrix(calibrates['alpha'], calibrates['beta'], calibrates['rho'], 
                                              calibrates['nu'], data['F'], data['K'], data['expiries'], method)
            results[fix+'_params']=calibrates['params']
    
    elif mark=='vega':
        # Test vega for over-specification analysis
        for fix in [0.2,0.4,0.6]:
            sabr=SABR.SABR_model(labels['label_ten'],labels['label_exp'],data['num_strikes'],
                                 labels['label_strikes'],data['strike_spreads'])
            calibrates=calibration(starting_guess,data['F'], data['K'], data['expiries'], data['MKT'],3, fix, method)
            results[fix+'_vols']=sabr.SABR_vol_matrix(calibrates['alpha'], calibrates['beta'], calibrates['rho'], 
                                              calibrates['nu'], data['F'], data['K'], data['expiries'], method)
            results[fix+'_params']=calibrates['params']
  
    return results

def objfunc(par,F,K,time,MKT,method='Hagan'):
    sum_sq_diff = 0
    if K[0]<=0:
        shift = 0.001 - K[0]
        for j in range(len(K)):
            K[j] = K[j] + shift
            F = F + shift
    for j in range(len(K)):
        if MKT[j] == 0:   
            diff = 0       
        elif F == K[j]:
            if method=='Hagan':
                V = (F*K[j])**((1-par[1])/2.)
                logFK = math.log(F/K[j])
                A = 1 + ( ((1-par[1])**2*par[0]**2)/(24.*(V**2)) + (par[0]*par[1]*par[3]*par[2])/(4.*V) + 
                         ((par[3]**2)*(2-3*(par[2]**2))/24.) ) * time
                B = 1 + (1/24.)*(((1-par[1])*logFK)**2) + (1/1920.)*(((1-par[1])*logFK)**4)
                VOL = (par[0]/V)*A
                diff = VOL - MKT[j]
            elif method=='Obloj':
                logFK = math.log(F/K[j])
                one_beta=1-par[1]
                one_betasqr=one_beta*one_beta
                fK=F*K[j]
                fK_beta=math.pow(fK,one_beta/2.0) 
                sigma_exp=(one_betasqr/24.0*par[0]*par[0]/fK_beta/fK_beta+0.25*par[2]*par[1]*par[3]*par[0]/fK_beta+
                           (2.0-3.0*par[2]*par[2])/24.0*par[3]*par[3])
                sigma=par[0]*math.pow(K[j],-one_beta)
                VOL=sigma*(1.0+sigma_exp*time)
                diff=VOL-MKT[j]
                
        elif F != K[j]: 
            if method=='Hagan':
                V = (F*K[j])**((1-par[1])/2.)
                logFK = math.log(F/K[j])
                z = (par[3]/par[0])*V*logFK
                x = math.log( ( math.sqrt(1-2*par[2]*z+z**2) + z - par[2] ) / (1-par[2]) )
                A = 1 + ( ((1-par[1])**2*par[0]**2)/(24.*(V**2)) + (par[0]*par[1]*par[3]*par[2])/(4.*V) + 
                         ((par[3]**2)*(2-3*(par[2]**2))/24.) ) * time
                B = 1 + (1/24.)*(((1-par[1])*logFK)**2) + (1/1920.)*(((1-par[1])*logFK)**4)
                VOL = (par[3]*logFK*A)/(x*B)
                diff=VOL-MKT[j]  
                    
            elif method=='Obloj': 
                logFK = math.log(F/K[j])
                one_beta=1-par[1]
                one_betasqr=one_beta*one_beta
                fK=F*K[j]
                fK_beta=math.pow(fK,one_beta/2.0) 
                sigma_exp=(one_betasqr/24.0*par[0]*par[0]/fK_beta/fK_beta+0.25*par[2]*par[1]*par[3]*par[0]/fK_beta+
                           (2.0-3.0*par[2]*par[2])/24.0*par[3]*par[3])
                if par[3]==0:
                    sigma=(1-par[1])*par[0]*logFK/(math.pow(F,(1-par[1]))-math.pow(K[j],(1-par[1])))
                elif par[1]==1:
                    z=par[3]*logFK/par[0]
                    sigma=par[3]*logFK/math.log((math.sqrt(1-2*par[2]*z+z*z)+z-par[2])/(1-par[2]))
                #elif (par[1]>0) and (par[1]<1):
                else:
                    z=par[3]*(math.pow(F,(1-par[1]))-math.pow(K[j],(1-par[1])))/par[0]/(1-par[1])
                    sigma=par[3]*logFK/math.log((math.sqrt(1-2*par[2]*z+z*z)+z-par[2])/(1-par[2]))
                VOL=sigma*(1.0+sigma_exp*time) 
                diff=VOL-MKT[j]
                   
        sum_sq_diff = sum_sq_diff + diff**2  
        obj = math.sqrt(sum_sq_diff)
    return obj

def calibration(starting_par,F,K,time,MKT,fix_par='auto',fix_no='auto',method='Hagan'):    
    starting_guess = np.array([0.001,0.5,0,0.001])
    
    if fix_par=='auto':
        pass
    else:
        starting_guess[fix_par]=fix_no
            
    alpha=len(F)*[starting_guess[0]]
    beta=len(F)*[starting_guess[1]]
    rho=len(F)*[starting_guess[2]]
    nu=len(F)*[starting_guess[3]]
    jacmat=len(F)*[starting_guess[3]]
        
    for i in range(len(F)):
        x0 = starting_par
        bnds = ((0.001,None) , (0,1) , (-0.999,0.999) , (0.001,None))
        if fix_par=='auto':
            res = minimize(objfunc,x0,(F[i],K[i],time[i],MKT[i],method),bounds=bnds,method='SLSQP') # for a constrained minimization of multivariate scalar functions
        else:
            res=minimize(objfunc,x0,(F[i],K[i],time[i],MKT[i],method),bounds=bnds,constraints={'type':'eq','fun':lambda par: par[fix_par]-fix_no},method='SLSQP') # with equality constraints added to calibrate another three with one parameter fixed
            
        alpha[i] = res.x[0]
        beta[i] = res.x[1]
        rho[i] = res.x[2]
        nu[i] = res.x[3]
        jacmat[i]=res.jac
        
    jacmat=pd.DataFrame(jacmat)
    params=pd.DataFrame(data=[alpha,beta,rho,nu,list(F),list(time)],index=['alpha','beta','rho','nu','F','time'])
    return {'alpha':alpha,'beta':beta,'rho':rho,'nu':nu,'params':params,'jacmat':jacmat}

