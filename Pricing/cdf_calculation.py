from scipy.stats import lognorm
import numpy as np
import pandas as pd
import math
#from matplotlib import pyplot as plt
from Pricing.black_pricing import BSPricer_SABR

def numerical_cdf(alpha,F,K,expiry,isCall=1,r=0,h=0.0001,vol_method='Hagan',vol_dist='lognormal'):
    [beta,rho,nu]=[0.5,0,0.001]
    #[alpha,beta,rho,nu] = [self.alpha,self.beta,self.rho,self.nu]
    bs = BSPricer_SABR(beta,rho,nu)
    if vol_dist=='lognormal':
        price = bs.BS_matrix(alpha,F,K,expiry,isCall,r,vol_method)
        price_plus = bs.BS_matrix(alpha,F,K+h,expiry,isCall,r,vol_method)
        cdf = (price_plus-price)/h
    #elif vol_dist=='normal':     
    return cdf

def digital_option_value(F,K,isCall):
    if isCall:
        if F>K:
            return 1.0
        else:
            return 0.0
    else:
        if F>K:
            return 0.0
        else:
            return 1.0



# In[ ]:



