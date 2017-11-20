import pandas as pd
import numpy as np
from scipy.stats import norm
import math
from Pricing.SABR import SABR_model

class BSPricer_SABR:
    def __init__(self,beta,rho,nu):
        self.beta = beta
        self.rho = rho
        self.nu = nu
        self.option_value = []
    
    def dPlusBlack(self,F_0,K,expiry,vol,r=0):
        D=math.exp(-r*expiry)
        d_Plus=(math.log(F_0/K)+(r+0.5*vol*vol)*expiry)/vol/math.sqrt(expiry)
        return d_Plus

    def dMinusBlack(self,F_0,K,expiry,vol,r=0):
        d_Minus=self.dPlusBlack(F_0,K,expiry,vol,r)-vol*math.sqrt(expiry)
        return d_Minus
    
    def black(self,F_0,K,expiry,vol,isCall,r=0):
    if expiry*vol == 0.0:
        if isCall:
            option_value=max(F_0*math.exp(r*expiry)-K,0.0)
        else:
            option_value=max(K-F_0*math.exp(r*expiry),0.0)
    else:
        d1=self.dPlusBlack(F_0,K,expiry,vol,r)
        d2=self.dMinusBlack(F_0,K,expiry,vol,r)
        if isCall:
            option_value=F_0*norm.cdf(d1)-y*norm.cdf(d2)
        else:
            option_value=y*norm.cdf(-d2)-F_0*norm.cdf(-d1)

    return option_value
    
    def price_lognorm_ivol(self,alpha,F_0,K,expiry,r=0,isCall=1,vol_method='Hagan'):
        sabr = SABR_model(self.beta,self.rho,self.nu)
        [beta,rho,nu] = [self.beta,self.rho,self.nu]
        D=math.exp(-r*expiry)
        if vol_method=='Hagan':
            vol = sabr.ivol_Hagan_lognorm(alpha,F_0,K,expiry)
        elif vol_method=='Obloj':
            vol = sabr.ivol_Obloj(alpha,F_0,K,expiry)
        if expiry*vol==0.0:
            if isCall:
                self.option_value=max(F_0/D-K,0.0)
            else:
                self.option_value=max(K-F_0/D,0.0)
        else:
            d1=self.dPlusBlack(F_0,K,expiry,vol,r)
            d2=self.dPlusBlack(F_0,K,expiry,vol,r)
            if isCall:
                option_value =F_0*norm.cdf(d1)-K*D*norm.cdf(d2)
            else:
                option_value =K*D*norm.cdf(-d2)-F_0*norm.cdf(-d1)
        return option_value
    
    def BS_vector(self,alpha,F_0,K,expiry,isCall,r,vol_method,i):
        sabr=SABR_model(self.beta,self.rho,self.nu)
        value = []
        for j in range(len(K)):
            if K[j]<=0:
                sabr.shift(F_0,K)
            V = self.price_lognorm_ivol(alpha,F_0,K[j],expiry,isCall,r,vol_method)
            value.append(V)
        return value

    def BS_matrix(self,alpha,F_0,K,expiry,isCall,r,vol_method): #F_T,expiry are vector, vol,K are a matrix
        option_value=[]
        for i in range(len(F_0)):
            V_vector = self.BS_vector(alpha,F_0[i],K[i],expiry[i],isCall,r,vol_method,i)
            option_value.append(V_vector)
            value_matrix = np.array(option_value) #columns=label_strikes)
        return value_matrix
    
    def find_ivol(self, option_price,F_T,K,expiry,r=0):
        sigma=0.20 #initial guess of sigma
        while sigma<1:
            black_implied=black(F_T,K,expiry,sigma,1,r)
            if option_price-black_implied<0.0001: #set precision of 0.0001
                return sigma
            print(sigma,option_price,black_implied)
            sigma+=0.01
        return "failture to find the right ivol!"
