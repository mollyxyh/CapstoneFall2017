import random
import math


def drawTwoRandomNumbers(rho):
    '''
    draw a pair of correlated random numbers
    @var rho: SABR rho
    '''

    z1 = random.gauss(0, 1)
    y1 = random.gauss(0, 1)

    term1 = z1 * rho
    term2 = (y1 * math.pow((1.0 - math.pow(rho, 2.0)), 0.5))
    x2 = term1 + term2

    return [z1, x2]

def simulateSABRMonteCarloEuler(no_of_sim, no_of_steps, expiry, F_0, alpha_0, beta, rho, nu):
    '''
    Monte Carlo SABR using Euler scheme.
    @var no_of_sim: Monte Carlo paths
    @var no_of_steps: discretization steps required to reach the option expiry date
    @var expiry: option expiry in years
    @var F_0: forward interest rate
    @var alpha_0: SABR alpha at t=0
    @var beta : SABR beta
    @var rho : SABR rho
    @var nu: SABR nu
    '''

    # Step length in years
    dt = float(expiry) / float(no_of_steps)
    dt_sqrt = math.sqrt(dt)
    no_of_sim_counter = 0
    simulated_forwards = []

    while no_of_sim_counter < no_of_sim:
        F_t = F_0
        alpha_t = alpha_0
        no_of_steps_counter = 1

        while no_of_steps_counter <= no_of_steps:
            # Zero absorbing boundary used for all beta except for beta=0 and beta=1
            if ((beta > 0 and beta < 1) and F_t <= 0):
                F_t = 0
                no_of_steps_counter = no_of_steps + 1
            else:
                # Generate two correlated random numbers
                rand = drawTwoRandomNumbers(rho)

                # Simulate the forward interest rate using Euler scheme.
                # Use the absolute for the diffusion to avoid numerical issues if the forward interest rate goes into negative
                dW_F = dt_sqrt * rand[0]
                F_b = math.pow(abs(F_t), beta)
                F_t = F_t + alpha_t * F_b * dW_F

                # Simulate the stochastic volatility using Euler scheme
                dW_a = dt_sqrt * rand[1]
                alpha_t = (alpha_t + nu * alpha_t * dW_a)

            no_of_steps_counter += 1

        # At the end of each path, we store the forward interest rate in a list
        simulated_forwards.append(F_t)
        no_of_sim_counter = no_of_sim_counter + 1

    return simulated_forwards

def simulateSABRMonteCarloMilstein(no_of_sim, no_of_steps, expiry, F_0, alpha_0, beta, rho, nu):
    '''
    Monte Carlo SABR using Milstein scheme.
    @var no_of_sim: Monte Carlo paths
    @var no_of_steps: discretization steps required to reach the option expiry date
    @var expiry: option expiry in years
    @var F_0: forward interest rate
    @var alpha_0: SABR alpha at t=0
    @var beta : SABR beta
    @var rho : SABR rho
    @var nu: SABR nu
    '''
    # Step length in years
    dt = float(expiry) / float(no_of_steps)
    dt_sqrt = math.sqrt(dt)
    no_of_sim_counter = 0
    simulated_forwards = []

    while no_of_sim_counter < no_of_sim:
        F_t = F_0
        alpha_t = alpha_0
        no_of_steps_counter = 1

        while no_of_steps_counter <= no_of_steps:
            # Zero absorbing boundary used for all beta except for beta=0 and beta=1
            if ((beta > 0 and beta < 1) and F_t <= 0):
                F_t = 0
                no_of_steps_counter = no_of_steps + 1

            else:
                # Generate two correlated random numbers
                rand = drawTwoRandomNumbers(rho)

                # Simulate the forward interest rate using Milstein scheme.
                # Use the absolute for the diffusion to avoid numerical issues if the forward interest rate goes into negative
                dW_F = dt_sqrt * rand[0]
                F_b = math.pow(abs(F_t), beta)
                exp_F = 2.0*beta-1.0
                F_t = (F_t + alpha_t * F_b * dW_F +
                       0.5 * beta * math.pow(alpha_t, 2.0) * math.pow(abs(F_t), exp_F) * (rand[0] * rand[0] - 1.0) * dt)

                # Simulate the stochastic volatility using Milstein scheme
                dW_a = dt_sqrt * rand[1]
                nu_sqr = math.pow(nu, 2.0)
                alpha_t = (alpha_t + nu * alpha_t * dW_a + 0.5 * nu_sqr * alpha_t * (rand[1] * rand[1] - 1.0) * dt)

            no_of_steps_counter += 1

        # At the end of each path, we store the forward interest rate in a list
        simulated_forwards.append(F_t)
        no_of_sim_counter = no_of_sim_counter + 1

    return simulated_forwards