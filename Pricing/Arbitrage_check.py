
import numpy as np
import pandas as pd
import xlrd
from matplotlib import pyplot as plt
from Pricing.black_pricing import black, dPlusBlack, dMinusBlack
from Pricing.pdf_calculation import numerical_pdf

def arbitrage_check(K,expiry,pdf):
    # K, expiry are list; pdf is matrix 
    for i in range(len(expiry)):
        for j in range(len(K)):
            if pdf[i][j] == 'No data':
                print("expiry = %s,K = %s: No data",expiry[i],K[j])
            elif pdf[i][j]<=0:
                print("expiry = %s,K = %s: Arbitrage Opportunity",expiry[i],K[j])

