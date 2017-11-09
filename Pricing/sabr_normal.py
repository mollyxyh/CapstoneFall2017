import math
from sympy import cosh,exp,sqrt,sinh,cos,sin
from scipy.stats import norm

T_grid=[0.25,0.5,0.75,1,2,5,10]
def NormalApprox(F_0,K,expiry,alpha_0,beta,rho,nu,T_grid):
    b0 = initial_vol(F_0,K,alpha_0,beta)
    B=I(K-F_0,beta,rho,nu)/b0
    S_min=0.0 # assume no shifting applied
    C=(2*I(S_min-F_0,beta,rho,nu)-I(K-F_0))/b0
    sumterm=0
    for i in range(T_grid.index(expiry)):
        sumterm+=term(i,rho,nu,expiry,T_grid,B,C,b0)
    return max(F_0-K,0)+math.pow(q(F_0-K,beta,rho,nu),0.25)/sqrt(2*math.pi)*b0*sumterm

# The initial normal volatility to use when approximating SABR with normal base
def initial_vol(F_0,K,alpha_0,beta):
    b0=alpha_0*(1-beta)*(K-F_0)/(math.pow(K,1-beta)-math.pow(F_0,1-beta))
    return b0

# Approximation for normal SABR__q(X)
def q(X,F_0,K,alpha_0,beta,rho):
    b0 = initial_vol(F_0,K,alpha_0,beta)
    qX = 1-2*rho*nu/beta0*x+math.pow((nu/b0)*X,2)
    return qX

def I(x,beta,rho,nu):
    return beta/nu*math.log((sqrt(q(x,beta,rho,nu)+rho+nu/beta*x))/(1+rho))

def k(i,rho,nu,expiry,T_grid,B):
    return -1.0/8*nu*nu+(math.log(Phi(T(i,expiry,T_grid),B,rho,nu))-math.log(Phi(T(i-1,expiry,T_grid),B,rho,nu)))/(T(i)-T(i-1))

def erf(x):
    return 2*norm.cdf(x*sqrt(2))-1

def ierf(x,y): # erf(x+y*j)
    term2=exp(-x*x)/(2*math.pi*x)*(1-cos(2*x*y)+sin(2*x*y)*j)
    term3=0
    for n in range(10):
        term3+=exp(-n*n/4.0)/(n*n+4*x*x)*(F(n,x,y)+G(n,x,y)*j)
    return erf(x)+term2+2.0/math.pi*exp(-x*x)*term3

def f(x,rho,nu):
    return 0.25*((1+rho)**2*exp(-2*nu*x)+(1-rho)**2*exp(2*nu*x)+2*(1-rho*rho))

def kappa(t,z,rho,nu): #first-order approximation
    return -1.0/8*nu*nu+3.0/16*nu*nu*(1-rho*rho)*(1.0/f(0,rho,nu)+1.0/f(z,rho,nu))

def Phi(t,z,rho,nu): #first-order approximation
    return exp(-1.0/8*nu*nu*t+3.0/16*nu*nu*(1-rho*rho)*(1.0/f(0,rho,nu)+1.0/f(z,rho,nu))*t)

def T(i,expiry,T_grid):
    if expiry!=10:
        T_grid=T_grid[:T_grid.index(expiry)+1]
    return T_grid[i]

def term(i,rho,nu,expiry,T_grid,B,C,b0):
    return exp(-(1.0/8*nu*nu+k(i,rho,nu,expiry,T_grid,B))*T(i-1,expiry,T_grid))*Phi(T(i-1,expiry,T_grid),I(0,beta,rho,nu)/b0,rho,nu)*J(i,rho,nu,expiry,T_grid,B,C,b0)

def F(n,x,y):
    return 2*x-2*x*cosh(n*y)*cos(2*x*y)+n*sinh(n*y)*sin(2*x*y)

def G(n,x,y):
    return 2*x*cosh(n*y)*sin(2*x*y)+n*sinh(n*y)*cos(2*x*y)

def J(i,rho,nu,expiry,T_grid,B,C,b0):
    kk1=k(i,rho,nu,expiry,T_grid,B)
    kk2=k(i-1,rho,nu,expiry,T_grid,B)
    mul1=sqrt(math.pi)/4.0/sqrt(-kk1)
    mul2=sqrt(math.pi)/4.0/sqrt(-kk2)
    tt1=T(i,expiry,T_grid)
    tt2=T(i-1,expiry,T_grid)
    B_b0=abs(float(B)/b0)
    C_b0=abs(float(C)/b0)
    [erf1B,erf1B_]=[ierf(B_b0/sqrt(tt1),sqrt(kk1*tt1)),ierf(-B_b0/sqrt(tt1),sqrt(kk1*tt1))]
    [erf2B,erf2B_]=[ierf(B_b0/sqrt(tt2),sqrt(kk2*tt2)),ierf(-B_b0/sqrt(tt2),sqrt(kk2*tt2))]
    [erf1C,erf1C_]=[ierf(C_b0/sqrt(tt1),sqrt(kk1*tt1)),ierf(-C_b0/sqrt(tt1),sqrt(kk1*tt1))]
    [erf2C,erf2C_]=[ierf(C_b0/sqrt(tt2),sqrt(kk2*tt2)),ierf(-C_b0/sqrt(tt2),sqrt(kk2*tt2))]
    term1=exp(2*B_b0*sqrt(-kk1))*(erf1B-1)+exp(-2*B_b0*sqrt(-kk1))*(erf1B_+1)
    term2=exp(2*B_b0*sqrt(-kk2))*(erf2B-1)+exp(-2*B_b0*sqrt(-kk2))*(erf2B_+1)
    term3=exp(2*C_b0*sqrt(-kk1))*(erf1C-1)+exp(-2*C_b0*sqrt(-kk1))*(erf1C_+1)
    term4=exp(2*C_b0*sqrt(-kk2))*(erf2C-1)+exp(-2*C_b0*sqrt(-kk2))*(erf2C_+1)
    return mul1*term1-mul2*term2-mul1*term3+mul2*term4
