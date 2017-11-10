import math
import cmath
from sympy import integrate
import sympy

def antonovLogNormalApprox(y, expiry, F_0, alpha_0, beta, nu, rho):
    '''
    function that returns the caplet Black price under Antonov approximation.
    @var y: option strike
    @var expiry: option expiry in years
    @var F_0: forward interest rate
    @var alpha_0: SABR alpha at t=0
    @var beta : SABR beta
    @var rho : SABR rho
    @var nu: SABR nu
    '''
    one_beta=1-beta
    q=cmath.pow(y,one_beta)/one_beta
    q0=cmath.pow(F_0,one_beta)/one_beta
    q_=(cmath.pow(y, one_beta)-cmath.pow(F_0,one_beta))/one_beta
    eta=1/abs(2.0*one_beta)

    nu_=cmath.sqrt(nu*nu-1.5*(nu*nu*rho*rho+alpha_0*nu*rho*one_beta/pow(F_0,one_beta)))
    p=phi(y,F_0,alpha_0,beta,nu,rho,nu_)
    alpha_0_=2*p*q_*nu_/(p*p-1)
    B_=B_min(y,F_0,alpha_0,beta,nu,rho,nu_) #complex no.
    alpha_1_=alpha_0_*nu_*nu_*(
        (0.5*cmath.log(alpha_0_*cmath.sqrt(q_*nu_*nu_+alpha_0_*alpha_0_))-B_)/cmath.log(p)/(p*p-1)*(
               p*p+1))
    alpha_=alpha_0_+expiry*alpha_1_

    s_minus=cmath.asinh(nu_*abs(q-q0)/alpha_)
    s_plus=cmath.asinh(nu_*abs(q+q0)/alpha_)

    term1=integrate(sympy.sin(eta*kappa(s,s_minus,s_plus))/sympy.sinh(s)*kernalG(expiry*nu_*nu_,s),
                      (s,s_minus,s_plus))
    term2=integrate(cmath.exp(-eta*psi(s,s_minus,s_plus))/cmath.sinh(s)*kernalG(expiry*nu_*nu_,s),
                      (s,s_plus,float("inf")))
    bone=term1+cmath.sin(eta*math.pi)*term2
    blk=max(F_0-y,0)+2.0/cmath.pi*cmath.sqrt(y*F_0)*bone
    return blk

def kappa(s,s_minus,s_plus):
    term1=sympy.sinh(s)*sympy.sinh(s)-sympy.sinh(s_minus)*sympy.sinh(s_minus)
    term2=sympy.sinh(s_plus)*sympy.sinh(s_plus)-sympy.sinh(s)*sympy.sinh(s)
    output=2.0*sympy.atan(sympy.sqrt(term1/term2))
    return output

def psi(s,s_minus,s_plus):
    term1=sympy.sinh(s)*sympy.sinh(s)-sympy.sinh(s_plus)*sympy.sinh(s_plus)
    term2=sympy.sinh(s)*sympy.sinh(s)-sympy.sinh(s_minus)*sympy.sinh(s_minus)
    output=2.0*sympy.atanh(sympy.sqrt(term1/term2))
    return output

def kernalG(tau,s):
    return sympy.sqrt(sympy.sinh(s)/s)*sympy.exp(-s*s/2.0/tau-tau/8.0)*(R(tau,s)+deltaR(tau,s))

def R(tau,s): # Approximation
    term1=1+3.0*tau*g(s)/8/s/s-5.0*tau*tau*(-8*s*s+3*pow(g(s),2)+24*g(s))/128/s**4
    term2=35.0*tau**3*(-40*s*s+3*pow(g(s),3)+24*pow(g(s),2)+120*g(s))/1024/s**6
    return term1+term2

def deltaR(tau,s):
    return sympy.exp(tau/8.0)-(3072+384*tau+24*tau*tau+tau**3)/3072.0
    
def g(s):
    return s*coth(s)-1

def coth(s):
    return sympy.cosh(s)/sympy.sinh(s)

def phi(y,F_0,alpha_0,beta,nu,rho,nu_):
    one_beta=1-beta
    q=pow(y,one_beta)/one_beta
    q_=(pow(y,one_beta)-pow(F_0,one_beta))/one_beta
    alpha_min=cmath.sqrt(nu*nu*q_*q+2.0*rho*nu*q_*alpha_0+alpha_0*alpha_0)
    base=(alpha_min+rho*alpha_0+nu*rho*q)/(1+rho)/alpha_0
    output=pow(base,nu_/nu)
    return output

def B_min(y,F_0,alpha_0,beta,nu,rho,nu_):
    one_beta=1-beta
    q=pow(y,one_beta)/one_beta
    q_=(pow(y,one_beta)-pow(F_0,one_beta))/one_beta
    alpha_min=cmath.sqrt(nu*nu*q_*q+2.0*rho*nu*q_*alpha_0+alpha_0*alpha_0)
    u0=(q_*nu*rho+alpha_0-alpha_min)/q_/nu/cmath.sqrt(1-rho*rho)
    L=alpha_min*one_beta/pow(y,one_beta)/nu/cmath.sqrt(1-rho*rho)
    if L.real>1:  # check for formula!!!
        I=1.0/cmath.sqrt(1-L*L)*cmath.log(
            (u0*(L+cmath.sqrt(L*L-1))+1)/(u0*(L-cmath.sqrt(L*L-1))+1)) 
    else:
        I=2.0/cmath.sqrt(1-L*L)*(
        cmath.atan((u0+L)/cmath.sqrt(1-L*L))-cmath.atan(L/cmath.sqrt(1-L*L)))
    phi0=cmath.acos(-(q_*nu+alpha_0*rho)/alpha_min)
    output=-0.5*beta/one_beta*rho/cmath.sqrt(1-rho*rho)*(cmath.pi-phi0-cmath.acos(rho-I)) #complex no.
    return output