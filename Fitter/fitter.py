import numpy as np
import pandas as pd
import math
#from sympy import solveset
#from sympy import Symbol, Eq
from scipy.optimize import minimize
from Pricing.SABR import SABR_model
from Pricing.Balland_sabr import Balland
from Pricing.Antonov_sabr import Antonov

class Fitter:
    def __init__(self, input_file):
        data = pd.read_excel('../Inputs/' + input_file)
        data.reset_index(inplace=True)
        K_spread = data.iloc[0, 3:].tolist()
        expiry = data.iloc[1:, 1].tolist()
        tenor = data.iloc[1:, 0].tolist()
        F = data.iloc[1:, 2].tolist()
        K = np.zeros((len(F), len(K_spread)))
        for i in range(len(F)):
            for j in range(len(K_spread)):
                K[i][j] = F[i] + 0.0001 * (K_spread[j])
        MKT = data.iloc[1:, 3:].values.tolist()
        self.MKT = MKT
        self.F = F
        self.K = K
        self.expiry = expiry
        self.tenor = tenor
        self.K_spread = K_spread
        
    def objfunc(self, par, F, K, expiry, MKT, method='Hagan_ln'):
        res = 0
        if K[0] <= 0: # shifting applied
            shift = 0.001 - K[0]
            for j in range(len(K)):
                K[j] = K[j] + shift
                F = F + shift
        if method == 'Hagan_ln':
            sabr = SABR_model(par[0],par[1], par[2])
            res_atm = minimize(self.objfunc_atm,0,(par[0],par[1], par[2], F, expiry, MKT[4], method))
            for j in range(len(K)):
                res += (sabr.ivol_Hagan_ln(res_atm.x[0],F,K[j],expiry) - MKT[j])**2
        elif method == 'Hagan_norm':
            sabr = SABR_model(par[0],par[1], par[2])
            res_atm = minimize(self.objfunc_atm,0,(par[0],par[1], par[2], F, expiry, MKT[4], method))
            for j in range(len(K)):
                res += (sabr.ivol_Hagan_norm(res_atm.x[0],F,K[j],expiry) - MKT[j])**2
        elif method == 'Obloj':
            sabr = SABR_model(par[0],par[1], par[2])
            res_atm = minimize(self.objfunc_atm,0,(par[0],par[1], par[2], F, expiry, MKT[4], method))
            for j in range(len(K)):
                res += (sabr.ivol_Obloj(res_atm.x[0],F,K[j],expiry) - MKT[j])**2
        #elif method == 'Balland':
            #sabr = Balland(par[1], par[2], par[3])
            #for j in range(len(K)):
                #res += (sabr.ivol(K[j],expiry,F,par[0],T_grid=[0.25,0.5,0.75,1,2,5,10]) - MKT[j])**2
        #elif method == 'Antonov':
            #sabr = Antonov(par[1], par[2], par[3])
            #for j in range(len(K)):
                #res += (sabr.ivol(K[j],expiry,F,par[0]) - MKT[j])**2
        obj = math.sqrt(res)
        return obj
    
    def objfunc_atm(self,alpha,beta,rho,nu,F,expiry, MKT, method='Hagan_ln'):
        sabr = SABR_model(beta,rho,nu)
        if method=='Hagan_ln':
            res=(sabr.ivol_Hagan_ln(alpha,F,F,expiry)-MKT)**2
        elif method=='Hagan_norm':
            res=(sabr.ivol_Hagan_norm(alpha,F,F,expiry)-MKT)**2
        elif method=='Obloj':
            res=(sabr.ivol_Obloj(alpha,F,F,expiry)-MKT)**2
        return res       

    def calibration(self, starting_par=np.array([0.5, 0, 0.001]), method='Hagan_ln', eqc='none'):
        [F, K, expiry, MKT, tenor] = [self.F, self.K, self.expiry, self.MKT, self.tenor]
        starting_guess = starting_par
        if eqc == 'none':
            pass
        else:
            starting_guess[eqc[0]] = eqc[1]
        alpha = len(F) * [0]
        beta = len(F) * [starting_guess[0]]
        rho = len(F) * [starting_guess[1]]
        nu = len(F) * [starting_guess[2]]
        jacmat = len(F) * [starting_guess[2]]

        for i in range(len(F)):
            x0 = starting_guess
            bnds = ((0, 1), (-0.999, 0.999), (0.001, None))
            if eqc == 'none':
                res = minimize(self.objfunc, x0, (F[i], K[i], expiry[i], MKT[i], method), bounds=bnds, method='SLSQP')
            elif len(eqc[0])==1:
                res = minimize(self.objfunc, x0, (F[i], K[i], expiry[i], MKT[i], method), bounds=bnds,
                               constraints={'type': 'eq', 'fun': lambda par: par[eqc[0][0]] - eqc[1]}, method='SLSQP')
            elif len(eqc[0])==2:
                res = minimize(self.objfunc, x0, (F[i], K[i], expiry[i], MKT[i], method), bounds=bnds,
                               constraints=[{'type': 'eq', 'fun': lambda par: par[eqc[0][0]] - eqc[1][0]},{'type': 'eq', 'fun': lambda par: par[eqc[0][1]] - eqc[1][1]}], method='SLSQP')

            #alpha[i] = res.x[0]
            beta[i] = res.x[0]
            rho[i] = res.x[1]
            nu[i] = res.x[2]
            jacmat[i] = res.jac
            
            res_atm = minimize(self.objfunc_atm,0,(beta[i],rho[i],nu[i], F[i], expiry[i], MKT[i][4], method))
            alpha[i] = res_atm.x[0]
    
            #ATM verification by scipy.optimize.minimize
            #res = minimize(self.objfunc_atm, alpha[i], (beta[i],rho[i],nu[i],F[i],expiry[i],MKT[i][4],method))
            #alpha[i] = res.x[0]
            
        jacmat = pd.DataFrame(jacmat)
        params = pd.DataFrame(data=[tenor,expiry,F,alpha,beta,rho,nu],index=['tenor','expiry','F','alpha','beta','rho','nu'])
        params = params.T
        return {'alpha': alpha, 'beta': beta, 'rho': rho, 'nu': nu, 'params': params, 'jacmat': jacmat}
    
    def ivol_SABR(self,alpha,beta,rho,nu,method='Hagan_ln'):
        sabr=SABR_model(0.5,0,0.25) #random nos
        [F,K,expiry]=[self.F,self.K,self.expiry]
        df=sabr.ivol_matrix(alpha,beta,rho,nu,F,K,expiry,method)
        df.columns=[-150,-100,-50,-25,0,25,50,100,150]
        df['F']=self.F
        df['tenor']=self.tenor
        df['expiry']=self.expiry
        df=df[['tenor','expiry','F',-150,-100,-50,-25,0,25,50,100,150]]
        return df