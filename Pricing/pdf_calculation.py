from scipy.stats import lognorm
from Pricing.black_pricing import black,dPlusBlack,dMinusBlack,bs_matrix

def numerical_pdf(bs_matrix):
    # price is a matrix
    m=[-150,-100,-50,-25,0,25,50,100,150]
    for i in len(F_T):
        for j in len(m):
            if j==0 or j==len(m)-1:
                num_pdf[i][j] = 'no data'
            else:
                num_pdf[i][j] = (bs_matrix[i][j+1]-2*bs_matrix[i][j]+bs_matrix[i][j-1])/(m[j]-m[j-1])**2
    return num_pdf

def lognormal_pdf(F,vol):
    lognormal_pdf = lognorm.pdf(F,vol)
    return lognormal_pdf
