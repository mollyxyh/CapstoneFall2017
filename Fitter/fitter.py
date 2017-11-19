import numpy as np
import pandas as pd
import math
from scipy.optimize import minimize
from Pricing.SABR import SABR_model

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

    def objfunc(self, par, F, K, expiry, MKT, method='Hagan'):
        sabr = SABR_model(par[1], par[2], par[3])
        res = 0
        if K[0] <= 0: # shifting applied
            shift = 0.001 - K[0]
            for j in range(len(K)):
                K[j] = K[j] + shift
                F = F + shift
        if method == 'Hagan':
            for j in range(len(K)):
                res += (sabr.ivol_Hagan(par[0],F,K[j],expiry) - MKT[j])**2
        elif method == 'Obloj':
            for j in range(len(K)):
                res += (sabr.ivol_Obloj(par[0],F,K[j],expiry) - MKT[j])**2
        obj = math.sqrt(res)
        return obj

    def calibration(self, starting_par=np.array([0.001, 0.5, 0, 0.001]), method='Hagan', eqc='none'):
        [F, K, expiry, MKT, tenor] = [self.F, self.K, self.expiry, self.MKT, self.tenor]
        starting_guess = starting_par
        if eqc == 'none':
            pass
        else:
            starting_guess[eqc[0]] = eqc[1]
        alpha = len(F) * [starting_guess[0]]
        beta = len(F) * [starting_guess[1]]
        rho = len(F) * [starting_guess[2]]
        nu = len(F) * [starting_guess[3]]
        jacmat = len(F) * [starting_guess[3]]

        for i in range(len(F)):
            x0 = starting_guess
            bnds = ((0.001, None), (0, 1), (-0.999, 0.999), (0.001, None))
            if eqc == 'none':
                res = minimize(self.objfunc, x0, (F[i], K[i], expiry[i], MKT[i], method), bounds=bnds, method='SLSQP')
            else:
                res = minimize(self.objfunc, x0, (F[i], K[i], expiry[i], MKT[i], method), bounds=bnds,
                               constraints={'type': 'eq', 'fun': lambda par: par[eqc[0]] - eqc[1]}, method='SLSQP')

            alpha[i] = res.x[0]
            beta[i] = res.x[1]
            rho[i] = res.x[2]
            nu[i] = res.x[3]
            jacmat[i] = res.jac

        jacmat = pd.DataFrame(jacmat)
        params = pd.DataFrame(data=[tenor,expiry,F,alpha,beta,rho,nu],index=['tenor','expiry','F','alpha','beta','rho','nu'])
        params = params.T
        return {'alpha': alpha, 'beta': beta, 'rho': rho, 'nu': nu, 'params': params, 'jacmat': jacmat}
    
    def ivol_SABR(self,alpha,beta,rho,nu,method='Hagan'):
        sabr=SABR_model(0.5,0,0.25) #random nos
        [F,K,expiry]=[self.F,self.K,self.expiry]
        df=sabr.ivol_matrix(alpha,beta,rho,nu,F,K,expiry,method)
        df.columns=[-150,-100,-50,-25,0,25,50,100,150]
        df['F']=self.F
        df['tenor']=self.tenor
        df['expiry']=self.expiry
        df=df[['tenor','expiry','F',-150,-100,-50,-25,0,25,50,100,150]]
        return df