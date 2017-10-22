
# coding: utf-8

# In[10]:

import numpy as np
import math
import csv
import xlrd
import matplotlib.pyplot as plt
with open('outprice.csv', newline='') as csvfile:
    pricereader = csv.reader(csvfile, delimiter=';',quotechar='|')
    pricelist = list(pricereader)
    print(pricelist)


# In[16]:

h = 0.0001*50
l1 = len(pricelist)

for k in range(2,l1):
    for j in range(2,len(pricelist[k])-1):
        pricelist[k][j] = float(pricelist[k][j])
        #print(pricelist[k][j])

sec_der_K = []
for k in range(2,l1):
    for j in range(0,len(pricelist[k])-1):
        if j==0 or j==1:
            der = pricelist[k][j]
        elif j==2 or j==len(pricelist[k])-2:
            der = 'No data'
        else:
            der = (pricelist[k][j+1] - 2*pricelist[k][j] + pricelist[k][j-1])/h**2
        sec_der_K.append(der)
#formated_der = [math.ceil(ele*10000)/10000 for ele in sec_der_K]

#print(sec_der_K)


# In[17]:

labels = ['Tenor','Expiry','-150','-100','-50','-25','ATM','25','50','100','150']
lst_parser = len(labels)
new_sec_der_K = [sec_der_K[k:k+lst_parser] for k in range(0, len(sec_der_K)-1, lst_parser)]
print(new_sec_der_K)


# In[18]:

with open("outpdf.csv", "w") as f:
    writer = csv.DictWriter(f, fieldnames=labels)
    writer.writeheader()
    writer = csv.writer(f)
    writer.writerows(new_sec_der_K)


# In[20]:

while True:
    try:
        file_input = xlrd.open_workbook('market_data.xlsx')     # load market data
    except:
        print('Input file is not in the directory!') 
    break
Market_data = file_input.sheet_by_name('Swaptions data')        # file input forward rates

strike_spreads=[]
j=0
while True:
    try:
        strike_spreads.append(int(Market_data.cell(1,3+j).value))
        j = j+1
    except:
        break
num_strikes = len(strike_spreads)

F = []
i=0
while True:
    try:
        F.append(Market_data.cell(2+i,2).value)
        i = i+1
    except:
        break
        
K = np.zeros((len(F),num_strikes))
for i in range(len(F)):
    for j in range(num_strikes):
        K[i][j] = F[i] + 0.0001*(strike_spreads[j])  
K_list = K.tolist()


# In[21]:

for k in range(len(K_list)):
    plt.plot(K_list[k][1:len(K_list[k])-1],new_sec_der_K[k][3:len(new_sec_der_K[k])-1],
             label= (new_sec_der_K[k][0]+new_sec_der_K[k][1]))
    plt.xlabel('strike price K')
    plt.ylabel('p.d.f of F_T')
    plt.show()


# In[ ]:



