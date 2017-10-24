import numpy as np
import pandas as pd
import math
from scipy.optimize import minimize
from Pricing.SABR import SABR_model

class Fitter:
    def __init__(self, input_file):
        self.input_file = input_file

    def input_read(self):
        data = pd.read_excel('../Inputs/' + self.input_file)
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

    def objfunc(self, par, F, K, expiry, MKT, method='Hagan'):
        res = 0
        if K[0] <= 0:
            shift = 0.001 - K[0]
            for j in range(len(K)):
                K[j] = K[j] + shift
                F = F + shift
        for j in range(len(K)):
            if MKT[j] == 0:  # no market data
                diff = 0
            elif F == K[j]:
                if method == 'Hagan':
                    V = (F * K[j]) ** ((1 - par[1]) / 2.)
                    logFK = math.log(F / K[j])
                    A = 1 + (((1 - par[1]) ** 2 * par[0] ** 2) / (24. * (V ** 2)) + (
                    par[0] * par[1] * par[3] * par[2]) / (4. * V) +
                             ((par[3] ** 2) * (2 - 3 * (par[2] ** 2)) / 24.)) * expiry
                    B = 1 + (1 / 24.) * (((1 - par[1]) * logFK) ** 2) + (1 / 1920.) * (((1 - par[1]) * logFK) ** 4)
                    ivol = (par[0] / V) * A
                    diff = ivol - MKT[j]
                elif method == 'Obloj':
                    logFK = math.log(F / K[j])
                    one_beta = 1 - par[1]
                    one_betasqr = one_beta * one_beta
                    fK = F * K[j]
                    fK_beta = math.pow(fK, one_beta / 2.0)
                    sigma_exp = (
                    one_betasqr / 24.0 * par[0] * par[0] / fK_beta / fK_beta + 0.25 * par[2] * par[1] * par[3] * par[
                        0] / fK_beta +
                    (2.0 - 3.0 * par[2] * par[2]) / 24.0 * par[3] * par[3])
                    sigma = par[0] * math.pow(K[j], -one_beta)
                    ivol = sigma * (1.0 + sigma_exp * expiry)
                    diff = ivol - MKT[j]
            elif F != K[j]:
                if method == 'Hagan':
                    V = (F * K[j]) ** ((1 - par[1]) / 2.)
                    logFK = math.log(F / K[j])
                    z = (par[3] / par[0]) * V * logFK
                    x = math.log((math.sqrt(1 - 2 * par[2] * z + z ** 2) + z - par[2]) / (1 - par[2]))
                    A = 1 + (((1 - par[1]) ** 2 * par[0] ** 2) / (24. * (V ** 2)) + (
                    par[0] * par[1] * par[3] * par[2]) / (4. * V) +
                             ((par[3] ** 2) * (2 - 3 * (par[2] ** 2)) / 24.)) * expiry
                    B = 1 + (1 / 24.) * (((1 - par[1]) * logFK) ** 2) + (1 / 1920.) * (((1 - par[1]) * logFK) ** 4)
                    ivol = (par[3] * logFK * A) / (x * B)
                    diff = ivol - MKT[j]
                elif method == 'Obloj':
                    logFK = math.log(F / K[j])
                    one_beta = 1 - par[1]
                    one_betasqr = one_beta * one_beta
                    fK = F * K[j]
                    fK_beta = math.pow(fK, one_beta / 2.0)
                    sigma_exp = (
                    one_betasqr / 24.0 * par[0] * par[0] / fK_beta / fK_beta + 0.25 * par[2] * par[1] * par[3] * par[
                        0] / fK_beta +
                    (2.0 - 3.0 * par[2] * par[2]) / 24.0 * par[3] * par[3])
                    if par[3] == 0:
                        sigma = (1 - par[1]) * par[0] * logFK / (
                        math.pow(F, (1 - par[1])) - math.pow(K[j], (1 - par[1])))
                    elif par[1] == 1:
                        z = par[3] * logFK / par[0]
                        sigma = par[3] * logFK / math.log(
                            (math.sqrt(1 - 2 * par[2] * z + z * z) + z - par[2]) / (1 - par[2]))
                    else:
                        z = par[3] * (math.pow(F, (1 - par[1])) - math.pow(K[j], (1 - par[1]))) / par[0] / (1 - par[1])
                        sigma = par[3] * logFK / math.log(
                            (math.sqrt(1 - 2 * par[2] * z + z * z) + z - par[2]) / (1 - par[2]))
                    ivol = sigma * (1.0 + sigma_exp * expiry)
                    diff = ivol - MKT[j]
            res += diff ** 2
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
        df.columns=['F',-150,-100,-50,-25,0,25,50,100,150]
        df['tenor']=self.tenor
        df['expiry']=self.expiry
        df=df[['tenor','expiry','F',-150,-100,-50,-25,0,25,50,100,150]]
        return df