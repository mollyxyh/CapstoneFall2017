import math
import numpy as np
import pandas as pd

class SABR_model:   
    def __init__(self,label_ten,label_exp,num_strikes,label_strikes,strike_spreads):
        self.label_ten=label_ten
        self.label_exp=label_exp
        self.num_strikes=num_strikes
        self.label_strikes=label_strikes
        self.strike_spreads=strike_spreads
        self.ivols=[]
        #self.ivols_diff=[]
        #self.parameters=[]
        
    def SABR(self,alpha,beta,rho,nu,F,K,time,method='Hagan'): # all variables are scalars
        if K <= 0:   # negative rates' problem, need to shift the smile
            VOL = 0
            #diff = 0
        elif F == K: # ATM formula
            if method=='Hagan':
                V = (F*K)**((1-beta)/2.)
                logFK = math.log(F/K)
                A = 1 + ( ((1-beta)**2*alpha**2)/(24.*(V**2)) + (alpha*beta*nu*rho)/(4.*V) + ((nu**2)*(2-3*(rho**2))/24.) ) * time
                B = 1 + (1/24.)*(((1-beta)*logFK)**2) + (1/1920.)*(((1-beta)*logFK)**4)
                VOL = (alpha/V)*A
                #diff = VOL - MKT                
            elif method=='Obloj':
                logFK = math.log(F/K)
                one_beta=1-beta
                one_betasqr=one_beta*one_beta
                fK=F*K
                fK_beta=math.pow(fK,one_beta/2.0) 
                sigma_exp=(one_betasqr/24.0*alpha*alpha/fK_beta/fK_beta+0.25*rho*beta*nu*alpha/fK_beta+(2.0-3.0*rho*rho)/24.0*nu*nu)
                sigma=alpha*math.pow(K,-one_beta)
                VOL=sigma*(1.0+sigma_exp*time)
                #diff=VOL-MKT      
        elif F != K: # not-ATM formula
            if method=='Hagan':
                V = (F*K)**((1-beta)/2.)
                logFK = math.log(F/K)
                z = (nu/alpha)*V*logFK
                x = math.log( ( math.sqrt(1-2*rho*z+z**2) + z - rho ) / (1-rho) )
                A = 1 + ( ((1-beta)**2*alpha**2)/(24.*(V**2)) + (alpha*beta*nu*rho)/(4.*V) + ((nu**2)*(2-3*(rho**2))/24.) ) * time
                B = 1 + (1/24.)*(((1-beta)*logFK)**2) + (1/1920.)*(((1-beta)*logFK)**4)
                VOL = (nu*logFK*A)/(x*B)
                #diff = VOL - MKT         
            elif method=='Obloj': ## Check for the formula!!!
                logFK = math.log(F/K)
                one_beta=1-beta
                one_betasqr=one_beta*one_beta
                fK=F*K
                fK_beta=math.pow(fK,one_beta/2.0) 
                sigma_exp=(one_betasqr/24.0*alpha*alpha/fK_beta/fK_beta+0.25*rho*beta*nu*alpha/fK_beta+(2.0-3.0*rho*rho)/24.0*nu*nu)
                if nu==0:
                    sigma=(1-beta)*alpha*logFK/(math.pow(F,(1-beta))-math.pow(K,(1-beta)))
                elif beta==1:
                    z=nu*logFK/alpha
                    sigma=nu*logFK/math.log((math.sqrt(1-2*rho*z+z*z)+z-rho)/(1-rho))
                else:
                #elif (beta>0) and (beta<1):
                    z=nu*(math.pow(F,(1-beta))-math.pow(K,(1-beta)))/alpha/(1-beta)
                    sigma=nu*logFK/math.log((math.sqrt(1-2*rho*z+z*z)+z-rho)/(1-rho))
                VOL=sigma*(1.0+sigma_exp*time) 
                #diff=VOL-MKT
        return VOL
          
    def smile(self,alpha,beta,rho,nu,F,K,time,i,method='Hagan'): # F, time and the parameters are scalars, K and MKT are vectors, i is the index for tenor/expiry label 
        temp=[self.label_ten[i],self.label_exp[i],F]
        #temp2=[self.label_ten[i],self.label_exp[i],F]
        for j in range(len(K)):
            if K[0]<= 0:
                self.shift(F,K)
            VOL = self.SABR(alpha,beta,rho,nu,F,K[j],time,method)
            temp.append(VOL)
            #temp2.append(diff)
        self.ivols.append(temp)
        #self.ivols_diff.append(temp2)
        #self.parameters.append([self.label_ten[i],self.label_exp[i],alpha,beta,rho,nu])
    
    def SABR_vol_matrix(self,alpha,beta,rho,nu,F,K,time,method='Hagan'): # F, time and the parameters are vectors, K and MKT are matrices
        label_strikes=['tenor','expiry','F']
        for j in range(len(self.label_strikes)):
            label_strikes.append(self.label_strikes[j])
        for i in range(len(F)):
            self.smile(alpha[i],beta[i],rho[i],nu[i],F[i],K[i],time[i],i,method)
        ivols=pd.DataFrame(data=self.ivols,columns=label_strikes) 
        #ivols_diff=pd.DataFrame(data=self.ivols_diff,columns=label_strikes) #,columns=label_strikes
        #params=pd.DataFrame(data=self.parameters,columns=['tenor','expiry','alpha','beta','rho','nu'])
        return ivols
            
    def shift(self,F,K):
        shift = 0.001 - K[0]
        for j in range(len(K)):
            K[j] = K[j] + shift
            F = F + shift
               

