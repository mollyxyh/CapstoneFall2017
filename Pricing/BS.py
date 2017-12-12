from scipy.stats import norm
from scipy.optimize import minimize
from Pricing.SABR import SABR_model
import math

def black(F_0,y,expiry,vol,isCall):
    if expiry*vol == 0.0:
        if isCall:
            option_value=max(F_0-y,0.0)
        else:
            option_value=max(y-F_0,0.0)
    else:
        d1=dPlusBlack(F_0=F_0,y=y,expiry=expiry,vol=vol)
        d2=dMinusBlack(F_0=F_0,y=y,expiry=expiry,vol=vol)
        if isCall:
            option_value=F_0*norm.cdf(d1)-y*norm.cdf(d2)
        else:
            option_value=y*norm.cdf(-d2)-F_0*norm.cdf(-d1)

    return option_value

def dPlusBlack(F_0,y,expiry,vol):
    d_plus=(math.log(F_0/y)+0.5*vol*vol*expiry)/vol/math.sqrt(expiry)
    return d_plus

def dMinusBlack(F_0,y,expiry,vol):
    d_minus=dPlusBlack(F_0=F_0,y=y,expiry=expiry,vol=vol)-vol*math.sqrt(expiry)
    return d_minus

def find_ivol(option_price,F_0,y,expiry):
        sigma=0.10 #initial guess of sigma
        while sigma<1:
            black_implied=black(F_0,y,expiry,sigma,1)
            if option_price-black_implied<0.000001: #set precision of 0.0001
                return sigma
            #print(sigma,option_price,black_implied)
            sigma+=0.0001
        return "failture to find the right ivol!"
    
def objfunc_atm(alpha,beta,rho,nu,F,expiry,MKT,method='Hagan_ln'):
    sabr = SABR_model(beta,rho,nu)
    if method=='Hagan_ln':
        res=(sabr.ivol_Hagan_ln(alpha,F,F,expiry)-MKT)**2
    elif method=='Hagan_norm':
        res=(sabr.ivol_Hagan_norm(alpha,F,F,expiry)-MKT)**2
    elif method=='Obloj':
        res=(sabr.ivol_Obloj(alpha,F,F,expiry)-MKT)**2
    return res