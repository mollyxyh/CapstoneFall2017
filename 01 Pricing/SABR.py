import xlrd
import math
import numpy as np
import pandas as pd
from scipy.optimize import minimize

class SABR_model:   
    def __init__(self,label_ten,label_exp,num_strikes,label_strikes,strike_spreads,outvol,vol_diff,parameters):
        self.label_ten=label_ten
        self.label_exp=label_exp
        self.num_strikes=num_strikes
        self.label_strikes=label_strikes
        self.strike_spreads=strike_spreads
        self.outvol=outvol
        self.vol_diff=vol_diff
        self.parameters=parameters
        
    def SABR(self,alpha,beta,rho,nu,F,K,time,MKT): # all variables are scalars
        if K <= 0:   # negative rates' problem, need to shift the smile
            VOL = 0
            diff = 0
        elif F == K: # ATM formula
            V = (F*K)**((1-beta)/2.)
            logFK = math.log(F/K)
            A = 1 + ( ((1-beta)**2*alpha**2)/(24.*(V**2)) + (alpha*beta*nu*rho)/(4.*V) + ((nu**2)*(2-3*(rho**2))/24.) ) * time
            B = 1 + (1/24.)*(((1-beta)*logFK)**2) + (1/1920.)*(((1-beta)*logFK)**4)
            VOL = (alpha/V)*A
            diff = VOL - MKT
        elif F != K: # not-ATM formula
            V = (F*K)**((1-beta)/2.)
            logFK = math.log(F/K)
            z = (nu/alpha)*V*logFK
            x = math.log( ( math.sqrt(1-2*rho*z+z**2) + z - rho ) / (1-rho) )
            A = 1 + ( ((1-beta)**2*alpha**2)/(24.*(V**2)) + (alpha*beta*nu*rho)/(4.*V) + ((nu**2)*(2-3*(rho**2))/24.) ) * time
            B = 1 + (1/24.)*(((1-beta)*logFK)**2) + (1/1920.)*(((1-beta)*logFK)**4)
            VOL = (nu*logFK*A)/(x*B)
            diff = VOL - MKT
        print round(VOL,4), '\t',
        self.outvol.write('%r;' %round(VOL,4) )
        if MKT==0:
            diff = 0
            self.vol_diff.write('%s;' %'No market data')
        else:
            self.vol_diff.write('%r;' %round(diff,4) )
            
    def smile(self,alpha,beta,rho,nu,F,K,time,MKT,i): # F, time and the parameters are scalars, K and MKT are vectors, i is the index for tenor/expiry label
        label_ten=self.label_ten
        label_exp=self.label_exp
        print label_ten[i] , '\t' , label_exp[i] , '\t' ,
        
        self.outvol.write('%s;%s;' %(label_ten[i],label_exp[i]))
        self.vol_diff.write('%s;%s;' %(label_ten[i],label_exp[i]))
        self.parameters.write('%s;%s;' %(label_ten[i],label_exp[i]))

        for j in range(len(K)):
            if K[0] <= 0:
                self.shift(F,K)
            self.SABR(alpha,beta,rho,nu,F,K[j],time,MKT[j])

        print ' '
        self.outvol.write('\n')
        self.vol_diff.write('\n')
        self.parameters.write('%f;%f;%f;%f;' %(alpha ,beta ,rho ,nu))
        self.parameters.write('\n')
    
    def SABR_vol_matrix(self,alpha,beta,rho,nu,F,K,time,MKT): # F, time and the parameters are vectors, K and MKT are matrices
        num_strikes=self.num_strikes
        label_strikes=self.label_strikes
        strike_spreads=self.strike_spreads
        
        print ' '
        print (2+((num_strikes-1)/2))*'       '+'SABR VOLATILITIES'
        print '  ' , '\t' , 'strikes:' , 
        for i in range(num_strikes):
            print label_strikes[i] , '\t' ,
        print ' '
        
        #outvol=self.outvol      # file output of volatilities
        #vol_diff=self.vol_diff  # file output differences between SABR and Market volatilities
        #parameters=self.parameters    # file output parameters
        
        self.outvol.write('%s;' %'SABR VOLATILITIES')
        self.outvol.write('\n')
        self.vol_diff.write('%s;' %'VOLATILITY DIFFERENCES')
        self.vol_diff.write('\n')
        self.parameters.write('%s;' %'PARAMETERS')
        self.parameters.write('\n')
        self.outvol.write('%s;%s;' %(' ','strikes:'))
        self.vol_diff.write('%s;%s;' %(' ','strikes:'))
        for j in range(len(strike_spreads)):
            self.outvol.write('%s;' %label_strikes[j])
            self.vol_diff.write('%s;' %label_strikes[j])
        self.outvol.write('\n')
        self.vol_diff.write('\n')
        print 'tenor' , '\t' ,   'expiry'
        self.parameters.write('%s;%s;%s;%s;%s;%s' %('tenor','expiry','alpha','beta','rho','nu'))
        self.parameters.write('\n')

        for i in range(len(F)):
            self.smile(alpha[i],beta[i],rho[i],nu[i],F[i],K[i],time[i],MKT[i],i)
            
    def shift(self,F,K):
        shift = 0.001 - K[0]
        for j in range(len(K)):
            K[j] = K[j] + shift
            F = F + shift
            
    def objfunc(self,par,F,K,time,MKT):
        sum_sq_diff = 0
        if K[0]<=0:
            self.shift(F,K)
        for j in range(len(K)):
            if MKT[j] == 0:   
                diff = 0       
            elif F == K[j]: 
                V = (F*K[j])**((1-par[1])/2.)
                logFK = math.log(F/K[j])
                A = 1 + ( ((1-par[1])**2*par[0]**2)/(24.*(V**2)) + (par[0]*par[1]*par[3]*par[2])/(4.*V) + ((par[3]**2)*(2-3*(par[2]**2))/24.) ) * time
                B = 1 + (1/24.)*(((1-par[1])*logFK)**2) + (1/1920.)*(((1-par[1])*logFK)**4)
                VOL = (par[0]/V)*A
                diff = VOL - MKT[j]
            elif F != K[j]: 
                V = (F*K[j])**((1-par[1])/2.)
                logFK = math.log(F/K[j])
                z = (par[3]/par[0])*V*logFK
                x = math.log( ( math.sqrt(1-2*par[2]*z+z**2) + z - par[2] ) / (1-par[2]) )
                A = 1 + ( ((1-par[1])**2*par[0]**2)/(24.*(V**2)) + (par[0]*par[1]*par[3]*par[2])/(4.*V) + ((par[3]**2)*(2-3*(par[2]**2))/24.) ) * time
                B = 1 + (1/24.)*(((1-par[1])*logFK)**2) + (1/1920.)*(((1-par[1])*logFK)**4)
                VOL = (par[3]*logFK*A)/(x*B)
                diff = VOL - MKT[j]  
            sum_sq_diff = sum_sq_diff + diff**2  
            obj = math.sqrt(sum_sq_diff)
        return obj
    
    def calibration(self,starting_par,F,K,time,MKT,fix_par='auto',fix_no='auto'):    
        global alpha,beta,rho,nu,jacmat
        starting_guess = np.array([0.001,0.5,0,0.001])  
        if fix_par=='auto':
            pass
        else:
            starting_guess[fix_par]=fix_no        
        alpha=len(F)*[starting_guess[0]]
        beta=len(F)*[starting_guess[1]]
        rho=len(F)*[starting_guess[2]]
        nu=len(F)*[starting_guess[3]]
        jacmat=len(F)*[starting_guess[3]]
        
        for i in range(len(F)):
            x0 = starting_par
            bnds = ((0.001,None) , (0,1) , (-0.999,0.999) , (0.001,None))
            if fix_par=='auto':
                res = minimize(self.objfunc,x0,(F[i],K[i],time[i],MKT[i]),bounds=bnds,method='SLSQP') # for a constrained minimization of multivariate scalar functions
            else:
                res = minimize(self.objfunc,x0,(F[i],K[i],time[i],MKT[i]),bounds=bnds,constraints={'type':'eq','fun':lambda par: par[fix_par]-fix_no},method='SLSQP') # for a constrained minimization of multivariate scalar functions
            
            alpha[i] = res.x[0]
            beta[i] = res.x[1]
            rho[i] = res.x[2]
            nu[i] = res.x[3]
            jacmat[i]=res.jac
            
        jacmat=pd.DataFrame(jacmat)
        params=pd.DataFrame(data=[alpha,beta,rho,nu],index=['alpha','beta','rho','nu'])
        if fix_par=='auto':
            jacmat.to_csv("../02 Fitter/parameters/jacmat.csv")     
            params.to_csv("../02 Fitter/parameters/params.csv")
        else:
            suffix='_'+str(fix_par)+'_'+str(fix_no)
            jacmat.to_csv("../02 Fitter/parameters/jacmat%s.csv"%(suffix))     
            params.to_csv("../02 Fitter/parameters/params%s.csv"%(suffix))
        return {'alpha':alpha,'beta':beta,'rho':rho,'nu':nu}
               