import math
from sympy import cosh, exp, sqrt, sinh, cos, sin
from scipy.stats import norm
from BS import *

T_grid = [0.25, 0.5, 0.75, 1, 2, 5, 10]

class Balland():
    def __init__(self,beta,rho,nu): # beta, rho and nu are stationary params of SABR model
        self.beta = beta
        self.rho = rho
        self.nu = nu
        
    def price(self, y, expiry, F_0, alpha_0, T_grid):
        b0 = alpha_0 * (1 - self.beta) * (y - F_0) / (math.pow(y, 1 - self.beta) - math.pow(F_0, 1 - self.beta))
        B = self.I(y - F_0) / b0
        S_min = 0 #assume no shifting applied
        C = (2 * self.I(S_min - F_0) - self.I(y - F_0)) / b0
        sumterm = 0
        for i in range(T_grid.index(expiry)):
            sumterm += self.term(i, expiry, T_grid, B, C, b0)
        return max(F_0 - y, 0) + math.pow(self.q(F_0 - y), 0.25) / sqrt(2 * math.pi) * b0 * sumterm
    
    def ivol(self, y, expiry, F_0, alpha_0, T_grid):
        p = self.price(y, expiry, F_0, alpha_0, T_grid)
        return find_ivol(p, F_0, y, expiry)
    
    def ivol_smile(self, alpha_0, F_0, y, expiry, i, T_grid): # alpha, F and expiry are scalars, K vectors, i the index for expiry
        ivols = []
        for j in range(len(y)):
            if y[0] <= 0:
                self.shift(F_0, y)
            ivol = self.ivol(y[j], expiry, F_0, alpha_0, T_grid)
            ivols.append(ivol)
        return ivols
    
    def ivol_matrix(self, alpha_0, beta, rho, nu, F_0, y, expiry, T_grid): # K is matrix, other variables are vectors
        ivols=[]
        for i in range(len(F_0)):
            self.set_params(beta[i], rho[i], nu[i])
            ivols.append(self.ivol_smile(alpha_0[i], F_0[i], y[i], expiry[i], i, T_grid))
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

    def erf(self, x):
        return 2 * norm.cdf(x * math.sqrt(2)) - 1

    def ierf(self, x, y):  # erf(x+y*j)
        term2 = exp(-x * x) / (2 * math.pi * x) * (1 - cos(2 * x * y) + sin(2 * x * y) * 1j)
        term3 = 0
        for n in range(10):
            term3 += exp(-n * n / 4.0) / (n * n + 4 * x * x) * (self.F(n, x, y) + self.G(n, x, y) * 1j)
        return self.erf(x) + term2 + 2.0 / math.pi * exp(-x * x) * term3

    def q(self, x):
        return 1 - 2 * self.rho * self.nu / self.beta * x + self.nu * self.nu / self.beta / self.beta * x * x

    def I(self, x):
        return self.beta / self.nu * math.log((sqrt(self.q(x) + self.rho + self.nu / self.beta * x)) / (1 + self.rho))

    def f(self, x):
        return 0.25 * ((1 + self.rho) ** 2 * exp(-2 * self.nu * x) + (1 - self.rho) ** 2 * exp(2 * self.nu * x) + 2 * (1 - self.rho * self.rho))

    def kappa(self, t, z):  # first-order approximation
        return -1.0 / 8 * self.nu * self.nu + 3.0 / 16 * self.nu * self.nu * (1 - self.rho * self.rho) * (1.0 / self.f(0) + 1.0 / self.f(z))

    def Phi(self, t, z):  # first-order approximation
        return exp(
            -1.0 / 8 * self.nu * self.nu * t + 3.0 / 16 * self.nu * self.nu * (1 - self.rho * self.rho) * (1.0 / self.f(0) + 1.0 / self.f(z)) * t)

    def T(self, i, expiry, T_grid):
        if expiry != 10:
            T_grid = T_grid[:T_grid.index(expiry) + 1]
        return T_grid[i]

    def k(self, i, expiry, T_grid, B):
        return -1.0 / 8 * self.nu * self.nu + (math.log(self.Phi(self.T(i, expiry, T_grid), B)) - math.log(
            self.Phi(self.T(i - 1, expiry, T_grid), B))) / (self.T(i, expiry, T_grid) - self.T(i - 1, expiry, T_grid))

    def term(self, i, expiry, T_grid, B, C, b0):
        return exp(-(1.0 / 8 * self.nu * self.nu + self.k(i, expiry, T_grid, B)) * self.T(i - 1, expiry, T_grid)) * self.Phi(
            self.T(i - 1, expiry, T_grid), self.I(0) / b0) * self.J(i, expiry, T_grid, B, C, b0)

    def F(self, n, x, y):
        return 2 * x - 2 * x * cosh(n * y) * cos(2 * x * y) + n * sinh(n * y) * sin(2 * x * y)

    def G(self, n, x, y):
        return 2 * x * cosh(n * y) * sin(2 * x * y) + n * sinh(n * y) * cos(2 * x * y)

    def J(self, i, expiry, T_grid, B, C, b0):
        kk1 = self.k(i, expiry, T_grid, B)
        kk2 = self.k(i - 1, expiry, T_grid, B)
        mul1 = sqrt(math.pi) / 4.0 / sqrt(-kk1)
        mul2 = sqrt(math.pi) / 4.0 / sqrt(-kk2)
        tt1 = self.T(i, expiry, T_grid)
        tt2 = self.T(i - 1, expiry, T_grid)
        B_b0 = abs(float(B) / b0)
        C_b0 = abs(float(C) / b0)
        [erf1B, erf1B_] = [self.ierf(B_b0 / math.sqrt(tt1), sqrt(kk1 * tt1)), self.ierf(-B_b0 / math.sqrt(tt1), sqrt(kk1 * tt1))]
        [erf2B, erf2B_] = [self.ierf(B_b0 / math.sqrt(tt2), sqrt(kk2 * tt2)), self.ierf(-B_b0 / math.sqrt(tt2), sqrt(kk2 * tt2))]
        [erf1C, erf1C_] = [self.ierf(C_b0 / math.sqrt(tt1), sqrt(kk1 * tt1)), self.ierf(-C_b0 / math.sqrt(tt1), sqrt(kk1 * tt1))]
        [erf2C, erf2C_] = [self.ierf(C_b0 / math.sqrt(tt2), sqrt(kk2 * tt2)), self.ierf(-C_b0 / math.sqrt(tt2), sqrt(kk2 * tt2))]
        term1 = exp(2 * B_b0 * sqrt(-kk1)) * (erf1B - 1) + exp(-2 * B_b0 * sqrt(-kk1)) * (erf1B_ + 1)
        term2 = exp(2 * B_b0 * sqrt(-kk2)) * (erf2B - 1) + exp(-2 * B_b0 * sqrt(-kk2)) * (erf2B_ + 1)
        term3 = exp(2 * C_b0 * sqrt(-kk1)) * (erf1C - 1) + exp(-2 * C_b0 * sqrt(-kk1)) * (erf1C_ + 1)
        term4 = exp(2 * C_b0 * sqrt(-kk2)) * (erf2C - 1) + exp(-2 * C_b0 * sqrt(-kk2)) * (erf2C_ + 1)
        return mul1 * term1 - mul2 * term2 - mul1 * term3 + mul2 * term4