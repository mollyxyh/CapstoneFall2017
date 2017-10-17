from scipy.stats import lognorm

def numerical_pdf(V,K,h):
    num_pdf = (V(K+h)-2*V(K)+V(K-h))/h**2
    return num_pdf

def lognormal_pdf(F,vol):
    lognormal_pdf = lognorm.pdf(F,vol)
    return lognormal_pdf
