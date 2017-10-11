import math


def oblojLogNormalApprox(y, expiry, F_0, alpha_0, beta, nu, rho):

    '''
    function that returns the Black implied volatility under Obloj lognormal approximation.
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

    sigma_exp = (one_betasqr / 24.0 * alpha_0 * alpha_0 / fK_beta / fK_beta + 0.25 * rho * beta * nu * alpha_0 / fK_beta + (2.0 - 3.0 * rho * rho) / 24.0 * nu * nu)

    if F_0 != y:
        if nu == 0:
            sigma1 = one_beta * alpha_0 * log_fK / (math.pow(F_0, one_beta) - math.pow(y, one_beta))
        elif beta == 1:
            z = nu * log_fK / alpha_0
            sigma1 = nu * log_fK / math.log((math.sqrt(1 - 2 * rho * z + z * z) + z - rho) / (1 - rho))
        elif (beta > 0) and (beta < 1):
            z = nu * (math.pow(F_0, one_beta) - math.pow(y, one_beta)) / alpha_0 / one_beta
            sigma1 = nu * log_fK / math.log((math.sqrt(1 - 2 * rho * z + z * z) + z - rho) / (1 - rho))
    else:
        sigma1 = alpha_0 * math.pow(y, -one_beta)

    sigma = sigma1 * (1.0 + sigma_exp * expiry)

    return sigma