from scipy.stats import norm
import math

def black(F_0,y,expiry,vol,isCall,r=0):
    option_value=0
    D=math.exp(-r*expiry)
    if expiry*vol==0.0:
        if isCall:
            option_value=max(F_0-y,0.0)
        else:
            option_value=max(y-F_0,0.0)
    else:
        d1=dPlusBlack(F_0,y,expiry,vol,r)
        d2=dPlusBlack(F_0,y,expiry,vol,r)
        if isCall:
            option_value=F_0*D*norm.cdf(d1)-y*D*norm.cdf(d2)
        else:
            option_value=y*D*norm.cdf(-d2)-F_0*D*norm.cdf(-d1)
    return option_value

def dPlusBlack(F_0,y,expiry,vol,r):
    d_Plus=(math.log(F_0/y*math.exp(-r*expiry))+(r+0.5*vol*vol)*expiry)/vol/math.sqrt(expiry)
    return d_Plus

def dMinusBlack(F_0,y,expiry,vol,r):
    d_minus=dPlusBlack(F_0,y,expiry,vol,r)-vol*math.sqrt(expiry)
    return d_minus

def find_ivol(option_price,F_0,y,expiry,r=0):
    sigma=0.20 #initial guess of sigma
    while sigma<1:
        black_implied=black(F_0,y,expiry,sigma,1,r)
        if option_price-black_implied<0.0001: #set precision of 0.0001
            return sigma
        print sigma,option_price,black_implied
        sigma+=0.01
    return "failture to find the right ivol!"
