import math
import numpy as np
import pandas as pd

class SABR_model:   
    def __init__(self,beta,rho,nu): # beta, rho and nu are stationary params of SABR model
        self.beta=beta
        self.rho=rho
        self.nu=nu
        
    def ivol_Hagan_lognorm(self,alpha,F,K,expiry): # alpha and F are state params of SABR model; all variables are scalars
        [beta,rho,nu]=[self.beta,self.rho,self.nu]
        if K<=0: # negative rates, shift needed
            ivol=0
        elif F==K: # ATM formula
            V=(F*K)**((1-beta)/2.)
            logFK=math.log(F/K)
            A=1+(((1-beta)**2*alpha**2)/(24.*(V**2))+(alpha*beta*nu*rho)/(4.*V)+((nu**2)*(2-3*(rho**2))/24.))*expiry
            B=1+(1/24.)*(((1-beta)*logFK)**2)+(1/1920.)*(((1-beta)*logFK)**4)
            ivol=(alpha/V)*A
        elif F!=K: # not-ATM formula
            V=(F*K)**((1-beta)/2.)
            logFK=math.log(F/K)
            z=(nu/alpha)*V*logFK
            x=math.log((math.sqrt(1-2*rho*z+z**2)+z-rho)/(1-rho))
            A=1+(((1-beta)**2*alpha**2)/(24.*(V**2))+(alpha*beta*nu*rho)/(4.*V)+((nu**2)*(2-3*(rho**2))/24.))*expiry
            B=1+(1/24.)*(((1-beta)*logFK)**2)+(1/1920.)*(((1-beta)*logFK)**4)
            ivol=(nu*logFK*A)/(x*B)       
        return ivol
    
    def zeta(self,alpha,F,K):
        [beta,rho,nu]=[self.beta,self.rho,self.nu]
        result = nu/alpha * (F*K)**((1-beta)/2)*math.log(F/K)
        return result
    
    def xhat(self,alpha,F,K):
        [beta,rho,nu]=[self.beta,self.rho,self.nu]
        zt = self.zeta(alpha,F,K)
        result = math.log((math.sqrt(1-2*rho*zt+zt**2)-rho+zt)/(1-rho))
        return result
        
    def ivol_Hagan_norm(self,alpha,F,K,expiry):
        [beta,rho,nu]=[self.beta,self.rho,self.nu]
        if K<=0: # negative rates, shift needed
            ivol=0
        elif F==K: # ATM formula (need to find the formula)
            V=(F*K)**((1-beta)/2.)
            logFK=math.log(F/K)
            A=1+(((1-beta)**2*alpha**2)/(24.*(V**2))+(alpha*beta*nu*rho)/(4.*V)+((nu**2)*(2-3*(rho**2))/24.))*expiry
            B=1+(1/24.)*(((1-beta)*logFK)**2)+(1/1920.)*(((1-beta)*logFK)**4)
            ivol=(alpha/V)*A
        elif F!=K: # not-ATM formula
            logFK=math.log(F/K)
            A=alpha*(F/K)**(beta/2.)
            B=(1+ 1/24.*(logFK**2) + 1/1920.*(logFK**4))/(1+ ((1-beta)**2)/24.*(logFK**2)+ ((1-beta)**4)/1920.*(logFK**4))
            C=self.zeta/self.xhat
            D=1+(((-beta*(2-beta)*alpha**2)/(24.*(FK)**(1-beta)) + (rho*alpha*nu*beta)/(4.*(FK)**((1-beta)/2.))+
                 ((2-3*rho**2)/24.)*(nu**2)))*expiry
            ivol = A*B*C*D     
        return ivol
    
    def ivol_Obloj(self,alpha,F,K,expiry): # alpha and F are state params of SABR model; all variables are scalars
        [beta,rho,nu]=[self.beta,self.rho,self.nu]
        if K<0 or K==0: # negative rates, shift needed
            ivol=0
        elif F==K: # ATM formula
            logFK=math.log(F/K)
            one_beta=1-beta
            one_betasqr=one_beta*one_beta
            fK=F*K
            fK_beta=math.pow(fK,one_beta/2.0) 
            sigma_exp=(one_betasqr/24.0*alpha*alpha/fK_beta/fK_beta+0.25*rho*beta*nu*alpha/fK_beta+(2.0-3.0*rho*rho)/24.0*nu*nu)
            sigma=alpha*math.pow(K,-one_beta)
            ivol=sigma*(1.0+sigma_exp*expiry)
        elif F!=K: # not-ATM formula
            logFK=math.log(F/K)
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
                z=nu*(math.pow(F,(1-beta))-math.pow(K,(1-beta)))/alpha/(1-beta)
                sigma=nu*logFK/math.log((math.sqrt(1-2*rho*z+z*z)+z-rho)/(1-rho))
            ivol=sigma*(1.0+sigma_exp*expiry) 
        return ivol
          
    def ivol_smile(self,alpha,F,K,expiry,i,method='Hagan'): # alpha, F and expiry are scalars, K vectors, i the index for expiry
        ivols=[]
        for j in range(len(K)):
            if K[0]<=0:
                self.shift(F,K)
            if method=='Hagan':
                ivol=self.ivol_Hagan_lognorm(alpha,F,K[j],expiry)
            elif method=='Obloj':
                ivol=self.ivol_Obloj(alpha,F,K[j],expiry)
            ivols.append(ivol)
        return ivols
    
    # TBD:
    def ivol_matrix_fit(self,alpha,beta,rho,nu,F,K,expiry,method='Hagan'): # K is matrix, other variables are vectors
        ivols=[]
        for i in range(len(F)):
            self.set_params(beta[i],rho[i],nu[i])
            ivols.append(self.ivol_smile(alpha[i],F[i],K[i],expiry[i],i,method))
        ivols=pd.DataFrame(data=ivols) 
        return ivols
    
    #New Function
    def ivol_matrix(self,alpha,F,K,expiry,method='Hagan'):
        ivols=[]
        for i in range(len(F)):
            ivols.append(self.ivol_smile(alpha,F[i],K[i],expiry[i],i,method))
        ivols=pd.DataFrame(data=ivols)
        return ivols
    
    def shift(self,F,K): # F and K are vectors
        shift=0.001-K[0]
        for j in range(len(K)):
            K[j]=K[j]+shift
            F=F+shift
            
    def set_params(self,beta,rho,nu):
        self.beta=beta
        self.rho=rho
        self.nu=nu
