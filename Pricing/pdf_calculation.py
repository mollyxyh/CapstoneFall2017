from scipy.stats import lognorm
import numpy as np
import pandas as pd
#import xlrd
import math
#from matplotlib import pyplot as plt
from Pricing.black_pricing import BSPricer_SABR
#from Pricing.Data_processor import data_reader,set_label,start_params

def numerical_pdf(alpha,beta,rho,nu,F,K,expiry,isCall=1,r=0,h=0.0001,vol_method='Hagan',vol_dist='lognormal'):
    #[beta,rho,nu]=[0.5,0,0.001]
    #[alpha,beta,rho,nu] = [self.alpha,self.beta,self.rho,self.nu]
    bs = BSPricer_SABR(beta,rho,nu)
    if vol_dist=='lognormal':
        price_minus = bs.BS_matrix(alpha,F,K-h,expiry,isCall,r,vol_method)
        price = bs.BS_matrix(alpha,F,K,expiry,isCall,r,vol_method)
        price_plus = bs.BS_matrix(alpha,F,K+h,expiry,isCall,r,vol_method)
        pdf = (price_plus-2*price+price_minus)/(h*h)
    #elif vol_dist=='normal':     
    return pdf

def lognormal_pdf(x,sigma,drift):
    a = 1/(x*sigma*math.sqrt(2*math.pi))
    b = math.exp(-(math.log(x)-drift)**2/(2*sigma*sigma))
    result = a * b 
    return result
