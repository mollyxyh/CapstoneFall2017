import math
import numpy as np
import sympy
from scipy.integrate import quad
from sympy.abc import s
from BS import *

class Antonov():
    def __init__(self,beta,rho,nu): # beta, rho and nu are stationary params of SABR model
        self.beta = beta
        self.rho = rho
        self.nu = nu
        
    def ivol(self, y, expiry, F_0, alpha_0):
        p = self.price(y, expiry, F_0, alpha_0)
        return find_ivol(p, F_0, y, expiry)
    
    def ivol_smile(self, alpha_0, F_0, y, expiry, i): # alpha, F and expiry are scalars, K vectors, i the index for expiry
        ivols = []
        for j in range(len(y)):
            if y[0] <= 0:
                self.shift(F_0, y)
            ivol = self.ivol(y[j], expiry, F_0, alpha_0)
            ivols.append(ivol)
        return ivols
    
    def ivol_matrix(self, alpha_0, beta, rho, nu, F_0, y, expiry, T_grid): # K is matrix, other variables are vectors
        ivols=[]
        for i in range(len(F_0)):
            self.set_params(beta[i], rho[i], nu[i])
            ivols.append(self.ivol_smile(alpha_0[i], F_0[i], y[i], expiry[i], i))
        ivols=pd.DataFrame(data=ivols) 
        return ivols
    
    def shift(self, F_0, y): # F and K are vectors
        shift=0.001-y[0]
        for j in range(len(y)):
            y[j] = y[j] + shift
            F_0 = F_0 + shift
            
    def set_params(self,beta,rho,nu):
        self.beta=beta
        self.rho=rho
        self.nu=nu

    def price(self, y, expiry, F_0, alpha_0):
        '''
        function that returns the caplet Black price under Antonov approximation.
        @var y: option strike
        @var expiry: option expiry in years
        @var F_0: forward interest rate
        @var alpha_0: SABR alpha at t=0
        '''
        [beta,rho,nu]=[self.beta,self.rho,self.nu]
        one_beta=1-beta
        q=math.pow(y,one_beta)/one_beta
        q0=math.pow(F_0,one_beta)/one_beta
        q_=(math.pow(y, one_beta)-math.pow(F_0,one_beta))/one_beta
        eta=1/abs(2.0*one_beta)

        nu_=math.sqrt(nu*nu-1.5*(nu*nu*rho*rho+alpha_0*nu*rho*one_beta/pow(F_0,one_beta)))
        p=self.phi(y,F_0,alpha_0)
        alpha_0_=2*p*q_*nu_/(p*p-1)
        B_=self.B_min(y,F_0,alpha_0)   
        R_=q_*nu_/alpha_0_
        alpha_min_=math.sqrt(nu_*nu_*q_*q_+alpha_0_*alpha_0_)
        alpha_min=math.sqrt(nu*nu*q_*q_+2.0*rho*nu*q_*alpha_0+alpha_0*alpha_0)
        alpha_1_=alpha_0_*nu_*nu_*math.sqrt(1+R_*R_)*(0.5*math.log(alpha_0*alpha_min/alpha_0_/alpha_min_)-B_)/R_/math.log(math.sqrt(1+R_*R_)+R_)
        alpha_=alpha_0_+expiry*alpha_1_

        s_minus=math.asinh(nu_*abs(q-q0)/alpha_)
        s_plus=math.asinh(nu_*abs(q+q0)/alpha_)

        pA_func1=sympy.lambdify([s],sympy.sin(eta*self.kappa(s,s_minus,s_plus))/sympy.sinh(s)*self.kernalG(expiry*nu_*nu_,s))
        term1=quad(pA_func1, s_minus, s_plus)[0]

        pA_func2=sympy.lambdify([s],sympy.exp(-eta*self.psi(s,s_minus,s_plus))/sympy.sinh(s)*self.kernalG(expiry*nu_*nu_,s))
        #term2=quad(pA_func2, s_plus, np.inf)[0] #math domain error rises
        term2=quad(pA_func2, s_plus, 20)[0] #replace np.inf with 20

        bone=term1+math.sin(eta*math.pi)*term2
        blk=max(F_0-y,0)+2.0/math.pi*math.sqrt(y*F_0)*bone
        return blk

    def kappa(self,s,s_minus,s_plus):
        term1=sympy.sinh(s)*sympy.sinh(s)-sympy.sinh(s_minus)*sympy.sinh(s_minus)
        term2=sympy.sinh(s_plus)*sympy.sinh(s_plus)-sympy.sinh(s)*sympy.sinh(s)
        output=2.0*sympy.atan(sympy.sqrt(term1/term2))
        return output

    def psi(self,s,s_minus,s_plus):
        term1=sympy.sinh(s)*sympy.sinh(s)-sympy.sinh(s_plus)*sympy.sinh(s_plus)
        term2=sympy.sinh(s)*sympy.sinh(s)-sympy.sinh(s_minus)*sympy.sinh(s_minus)
        output=2.0*sympy.atanh(sympy.sqrt(term1/term2))
        return output

    def kernalG(self,tau,s):
        return sympy.sqrt(sympy.sinh(s)/s)*sympy.exp(-s*s/2.0/tau-tau/8.0)*(self.R(tau,s)+self.deltaR(tau,s))

    def R(self,tau,s): # Approximation
        term1=1+3.0*tau*self.g(s)/8/s/s-5.0*tau*tau*(-8*s*s+3*pow(self.g(s),2)+24*self.g(s))/128/s**4
        term2=35.0*tau**3*(-40*s*s+3*pow(self.g(s),3)+24*pow(self.g(s),2)+120*self.g(s))/1024/s**6
        return term1+term2

    def deltaR(self,tau,s):
        return sympy.exp(tau/8.0)-(3072+384*tau+24*tau*tau+tau**3)/3072.0

    def g(self,s):
        return s*self.coth(s)-1

    def coth(self,s):
        return sympy.cosh(s)/sympy.sinh(s)

    def phi(self,y,F_0,alpha_0):
        [beta,rho,nu]=[self.beta,self.rho,self.nu]
        one_beta=1-beta
        q=pow(y,one_beta)/one_beta
        q_=(pow(y,one_beta)-pow(F_0,one_beta))/one_beta
        nu_=math.sqrt(nu*nu-1.5*(nu*nu*rho*rho+alpha_0*nu*rho*one_beta/pow(F_0,one_beta)))
        alpha_min=math.sqrt(nu*nu*q_*q+2.0*rho*nu*q_*alpha_0+alpha_0*alpha_0)
        base=(alpha_min+rho*alpha_0+nu*rho*q)/(1+rho)/alpha_0
        output=pow(base,nu_/nu)
        return output

    def B_min(self,y,F_0,alpha_0):
        [beta,rho,nu]=[self.beta,self.rho,self.nu]
        one_beta=1-beta
        nu_=math.sqrt(nu*nu-1.5*(nu*nu*rho*rho+alpha_0*nu*rho*one_beta/pow(F_0,one_beta)))
        q=pow(y,one_beta)/one_beta
        q_=(pow(y,one_beta)-pow(F_0,one_beta))/one_beta
        alpha_min=math.sqrt(nu*nu*q_*q+2.0*rho*nu*q_*alpha_0+alpha_0*alpha_0)
        u0=(q_*nu*rho+alpha_0-alpha_min)/q_/nu/math.sqrt(1-rho*rho)
        L=alpha_min*one_beta/pow(y,one_beta)/nu/math.sqrt(1-rho*rho)
        if L>1: 
            I=1.0/math.sqrt(L*L-1)*math.log(
                (u0*(L+math.sqrt(L*L-1))+1)/(u0*(L-math.sqrt(L*L-1))+1)) 
        elif L<1:
            I=2.0/math.sqrt(1-L*L)*(
            math.atan((u0+L)/math.sqrt(1-L*L))-math.atan(L/math.sqrt(1-L*L)))
        phi0=math.acos(-(q_*nu+alpha_0*rho)/alpha_min)
        output=-0.5*beta/one_beta*rho/math.sqrt(1-rho*rho)*(math.pi-phi0-math.acos(rho)-I)
        return output