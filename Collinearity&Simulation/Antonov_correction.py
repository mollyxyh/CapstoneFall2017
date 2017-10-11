import math
from sympy import integrate,cosh,exp,sqrt
from sympy.abc import s,x


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
    one_beta = 1 - beta
    one_betasqr = one_beta * one_beta
    fK = F_0 * y
    fK_beta = math.pow(fK, one_beta / 2.0)
    log_fK = math.log(F_0 / y)

    q = math.pow(y, one_beta) / one_beta
    q0 = math.pow(F_0, one_beta) / one_beta
    q_ = (math.pow(y, one_beta) - math.pow(F_0, one_beta)) / one_beta
    eta = 1 / abs(2.0 * one_beta)

    nu_ = math.sqrt(nu * nu - 3 / 2.0 * (nu * nu * rho * rho + alpha_0 * nu * rho * one_beta / math.pow(F_0, one_beta)))
    p = phi(y, F_0, alpha_0, beta, nu, rho, nu_)
    alpha_0_ = 2 * p * q_ * nu_ / (p * p - 1)
    B_ = B_min(y, F_0, alpha_0, beta, nu, rho, nu_)
    alpha_1_ = alpha_0_ * nu_ * nu_ * (-1.0) * (
    0.5 * math.log(alpha_0_ * math.sqrt(q_ * nu_ * nu_ + alpha_0_ * alpha_0_)) - B_) / math.log(p) / (p * p - 1) * (
               p * p + 1)
    alpha_ = alpha_0_ + expiry * alpha_1_

    s_minus = math.asinh(nu_ * abs(q - q0) / alpha_)
    s_plus = math.asinh(nu_ * abs(q + q0) / alpha_)

    term1 = integrate(math.sin(eta * kappa(s, s_minus, s_plus)) / math.sinh(s) * kernalG(expiry * nu_ * nu_, s),
                      (s, s_minus, s_plus))
    term2 = integrate(math.log(-eta * psi(s, s_minus, s_plus)) / math.sinh(s) * kernalG(expiry * nu_ * nu_, s),
                      (s, s_plus, float("inf")))
    bone = term1 + math.sin(eta * math.pi) * term2
    blk = max(F_0 - y, 0) + 2.0 / math.pi * math.sqrt(y * F_0) * bone
    return blk



def kappa(s,s_minus,s_plus):
    term1=math.asinh(s)*math.asinh(s)-math.asinh(s_minus)*math.asinh(s_minus)
    term2=math.asinh(s_plus)*math.asinh(s_plus)-math.asinh(s)*math.asinh(s)
    output=2.0*math.atan(math.sqrt(term1/term2))
    return output

def psi(s,s_minus,s_plus):
    term1=math.asinh(s)*math.asinh(s)-math.asinh(s_plus)*math.asinh(s_plus)
    term2=math.asinh(s)*math.asinh(s)-math.asinh(s_minus)*math.asinh(s_minus)
    output=2.0*math.atanh(math.sqrt(term1/term2))
    return output

def kernalG(tau,s):
    integral=integrate(x*exp(-x*x/tau)*sqrt(cosh(x)-cosh(s)),(x,s,float("inf")))
    return integral

def phi(y,F_0,alpha_0,beta,nu,rho,nu_):
    one_beta=1-beta
    q=math.pow(y,one_beta)/one_beta
    q_=(math.pow(y,one_beta)-math.pow(F_0,one_beta))/one_beta
    alpha_min=math.sqrt(nu*nu*q_*q+2.0*rho*nu*q_*alpha_0+alpha_0*alpha_0)
    base=(alpha_min+rho*alpha_0+nu*rho*q)/(1+rho)/alpha_0
    output=math.pow(base,nu_/nu)
    return output

def B_min(y, F_0, alpha_0, beta, nu, rho, nu_):
    one_beta = 1 - beta
    q = math.pow(y, one_beta) / one_beta
    q_ = (math.pow(y, one_beta) - math.pow(F_0, one_beta)) / one_beta
    alpha_min = math.sqrt(nu * nu * q_ * q + 2.0 * rho * nu * q_ * alpha_0 + alpha_0 * alpha_0)
    u0 = (q_ * nu * rho + alpha_0 - alpha_min) / q_ / nu / math.sqrt(1 - rho * rho)
    L = alpha_min * one_beta / math.pow(y, one_beta) / nu / math.sqrt(1 - rho * rho)
    if L > 1:  # check for formula!
        I = 1.0 / math.sqrt(L * L - 1) * math.log(
            (u0 * (L + math.sqrt(L * L - 1)) + 1) / (u0 * (L - math.sqrt(L * L - 1)) + 1))
    else:
        I = 2.0 / math.sqrt(L * L - 1) * (
        math.atan((u0 + L) / math.sqrt(1 - L * L)) - math.atan(L / math.sqrt(1 - L * L)))

    phi0 = math.acos(-(q_ * nu + alpha_0 * rho) / alpha_min)
    output = -0.5 * beta / one_beta * rho / math.sqrt(1 - rho * rho) * (math.pi - phi0 - math.acos(rho - I))
    return output