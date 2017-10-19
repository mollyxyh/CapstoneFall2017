from scipy.stats import norm
import math

def black(F_T,K,expiry,vol,isCall,r=0):
    #option_value=0
    D=math.exp(-r*expiry)
    if expiry*vol==0.0:
        if isCall:
            option_value=max(F_T-K,0.0)
        else:
            option_value=max(K-F_T,0.0)
    else:
        d1=dPlusBlack(F_T,K,expiry,vol,r)
        d2=dPlusBlack(F_T,K,expiry,vol,r)
        if isCall:
            option_value=F_T*D*norm.cdf(d1)-K*D*norm.cdf(d2)
        else:
            option_value=K*D*norm.cdf(-d2)-F_T*D*norm.cdf(-d1)
    return option_value

def dPlusBlack(F_T,K,expiry,vol,r):
    D=math.exp(-r*expiry)
    d_Plus=(math.log(F_T*D/K)+(r+0.5*vol*vol)*expiry)/vol/math.sqrt(expiry)
    return d_Plus

def dMinusBlack(F_T,K,expiry,vol,r):
    d_Minus=dPlusBlack(F_T,K,expiry,vol,r)-vol*math.sqrt(expiry)
    return d_Minus

def find_ivol(option_price,F_T,K,expiry,r=0):
    sigma=0.20 #initial guess of sigma
    while sigma<1:
        black_implied=black(F_T,K,expiry,sigma,1,r)
        if option_price-black_implied<0.0001: #set precision of 0.0001
            return sigma
        print sigma,option_price,black_implied
        sigma+=0.01
    return "failture to find the right ivol!"

def bs_vector(F_T,K,expiry,vol,isCall,i,r=0):
    for j in range(len(K)):
        black(F_T,K[j],expiry,vol[j],isCall,r=0)

def bs_matrix(F_T,K,expiry,vol,isCall,r=0): #F_T,expiry are vector, vol,K are a matrix
    for i in range(len(F_T)):
        bs_vector(F_T[i],K[i],expiry[i],vol[i],isCall,i,r=0)

