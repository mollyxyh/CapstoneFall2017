import math
import numpy as np
import pandas as pd
from scipy.integrate import quad,ode
import sympy
from sympy.abc import s

class ZABR_model:   
    def __init__(self,beta,rho,nu,gamma): #beta, rho and nu are stationary params of SABR model
        self.beta=beta
        self.rho=rho
        self.nu=nu
        self.gamma=gamma
        
    def ivols_ln(self,alpha,F,expiry,K): #K is vector
        [beta,rho,nu,gamma]=[self.beta,self.rho,self.nu,self.gamma]
        y=[]
        for i in range(len(K)):
            y.append((pow(F,1-beta)-pow(K[i],1-beta))/(1-beta)*pow(alpha,(gamma-2))) #analytical integral value
            #func=sympy.lambdify([s],1./self.vol(alpha,F))
            #y.append(quad(func,y,F)[0]*pow(alpha,gamma-2)) #approximated integral value
        
        f0,y0=0,0
        fy=[]
        for i in range(len(K)):
            solver=ode(self.odefunc).set_integrator('dopri5') #same algo as ode45 in MATLAB
            solver.set_initial_value(f0,y0)
            solver._integrator.iwork[2]=-1
            fy.append(solver.integrate(y[i],step=True))
            f0=fy[i]
            y0=y[i]
            
        x=[pow(alpha,1-gamma)*item for item in fy]
        ivol_ln=[]
        for i in range(len(K)):
            if F!=K[i]:
                ivol_ln.append(math.log(F/K[i])/x[i])
            else:
                ivol_ln.append(alpha/pow(F,1-beta))
        return ivol_ln
    
    def ivols_n(self,alpha,F,expiry,K): #K is vector
        [beta,rho,nu,gamma]=[self.beta,self.rho,self.nu,self.gamma]
        y=[]
        for i in range(len(K)):
            y.append((pow(F,1-beta)-pow(K[i],1-beta))/(1-beta)*pow(alpha,(gamma-2))) #analytical integral value
            #func=sympy.lambdify([s],1./self.vol(alpha,F))
            #y.append(quad(func,y,F)[0]*pow(alpha,gamma-2)) #approximated integral value
        
        f0,y0=0,0
        fy=[]
        for i in range(len(K)):
            solver=ode(self.odefunc).set_integrator('dopri5') #same algo as ode45 in MATLAB
            solver.set_initial_value(f0,y0)
            solver._integrator.iwork[2]=-1
            fy.append(solver.integrate(y[i],step=True))
            f0=fy[i]
            y0=y[i]
            
        x=[pow(alpha,1-gamma)*item for item in fy]
        ivol_n=[]
        for i in range(len(K)):
            if F!=K[i]:
                ivol_n.append((F-K[i])/x[i])
            else:
                ivol_n.append(pow(F,beta)*alpha)
        return ivol_n
        
    
    def A(self,y):
        [rho,nu,gamma]=[self.rho,self.nu,self.gamma]
        return 1+(gamma-2)**2*nu*nu*y*y+2*rho*(gamma-2)*nu*y
    
    def B(self,y):
        [rho,nu,gamma]=[self.rho,self.nu,self.gamma]
        return 2*rho*(1-gamma)*nu+2*(1-gamma)*(gamma-2)*nu*nu*y
    
    def odefunc(self,t,f):
        [nu,gamma]=[self.nu,self.gamma]
        C=(1-gamma)**2*nu*nu
        return 0.5*(-self.B(t)*f+math.sqrt(self.B(t)**2.*pow(f,2)-4*self.A(t)*(C*f*f-1)))/self.A(t)
    
    def vol(self,alpha,F):
        return alpha*pow(F,self.beta)

