import pandas as pd
import numpy as np
from scipy.stats import norm
import math
import Pricing.SABR as SABR
from Pricing.Data_processor import data_reader,set_label

class BS_pricing:   
    def __init__(self,label_ten,label_exp,num_strikes,label_strikes,strike_spreads):
        self.label_ten=label_ten
        self.label_exp=label_exp
        self.num_strikes=num_strikes
        self.label_strikes=label_strikes
        self.strike_spreads=strike_spreads
        self.option_value=[]

    def black(self,F_T,K,expiry,vol,isCall,r=0):
    #option_value=0
        D=math.exp(-r*expiry)
        if expiry*vol==0.0:
            if isCall:
                V=max(F_T-K,0.0)
            else:
                V=max(K-F_T,0.0)
        else:
            d1=self.dPlusBlack(F_T,K,expiry,vol,r)
            d2=self.dPlusBlack(F_T,K,expiry,vol,r)
            if isCall:
                V =F_T*D*norm.cdf(d1)-K*D*norm.cdf(d2)
            else:
                V =K*D*norm.cdf(-d2)-F_T*D*norm.cdf(-d1)
        return V

<<<<<<< HEAD
    def dPlusBlack(self,F_T,K,expiry,vol,r):
        D=math.exp(-r*expiry)
        d_Plus=(math.log(F_T*D/K)+(r+0.5*vol*vol)*expiry)/vol/math.sqrt(expiry)
        return d_Plus
=======
def find_ivol(option_price,F_T,K,expiry,r=0):
    sigma=0.20 #initial guess of sigma
    while sigma<1:
        black_implied=black(F_T,K,expiry,sigma,1,r)
        if option_price-black_implied<0.000001: #set precision of 0.000001
            return sigma
        print sigma,option_price,black_implied
        sigma+=0.01
    return "failture to find the right ivol!"
>>>>>>> 72a9f59ae4ee4b5c56bb35e2a4e6c9d075ae2716

    def dMinusBlack(self,F_T,K,expiry,vol,r):
        d_Minus=dPlusBlack(F_T,K,expiry,vol,r)-vol*math.sqrt(expiry)
        return d_Minus

    def bs_vector(self,F_T,K,expiry,vol,isCall,i,r=0):
        sabr=SABR.SABR_model(self.label_ten,self.label_exp,self.num_strikes,self.label_strikes,self.strike_spreads)
        temp=[self.label_ten[i],self.label_exp[i],F_T]
        for j in range(len(K)):
            if K[j]<=0:
                sabr.shift(F_T,K)
            V = self.black(F_T,K[j],expiry,vol.values[j],isCall,r=0)
            temp.append(V)
        self.option_value.append(temp)                                    

    def bs_matrix(self,F_T,K,expiry,vol,isCall,r=0): #F_T,expiry are vector, vol,K are a matrix
        label_strikes=['tenor','expiry','F']
        for j in range(len(self.label_strikes)):
            label_strikes.append(self.label_strikes[j])
        for i in range(len(F_T)):
            v = vol.loc[i][3:]
            option_value = self.bs_vector(F_T[i],K[i],expiry[i],v,isCall,i,r=0)
            option_value = pd.DataFrame(data=self.option_value,columns=label_strikes)
        return option_value
    
    def find_ivol(self,option_price,F_T,K,expiry,r=0):
        sigma=0.20 #initial guess of sigma
        while sigma<1:
            black_implied=black(F_T,K,expiry,sigma,1,r)
            if option_price-black_implied<0.0001: #set precision of 0.0001
                return sigma
            print(sigma,option_price,black_implied)
            sigma+=0.01
        return "failture to find the right ivol!"
