import numpy as np
import pandas as pd
import math
import numpy.matlib

# =============================================================================
#Make sure to have 'PAH_6MWD2.xlsx' and 'PAH_with_6MWD.xlsx' in the same direcctory
# =============================================================================
healthy_val = pd.read_excel('PAH_with_6MWD.xlsx')['Predicted 6MWD']
xls = pd.ExcelFile('PAH_6MWD2.xlsx')
sheet_names = xls.sheet_names

initial =[]
proj = []
health_mean = []
slope_fitbit = []
sample = []
sample_smooth = []
for i in range(len(sheet_names)): 


    df = pd.read_excel(xls, sheet_names[i])
    raw_data = df.iloc[:,2].dropna()
    smoothed_data =  raw_data.ewm(alpha=0.3).mean() ## Smoothed data
    m1, b1 = np.polyfit(range(1,len(smoothed_data)+1),smoothed_data, 1)
    
    data0 =smoothed_data.iloc[0]
    projection = m1*len(smoothed_data)+data0  ## calcultinng the projection of 6mwd based on the smoothed slope
    First_val = df.iloc[0,1]
    Second_val = df.iloc[-1,1]
        
    if  not (math.isnan(Second_val) or math.isnan(First_val)):
        
         initial.append(data0)
         proj.append(projection)  ### colecting the projected daat
         
         spli = sheet_names[i].split('PAH')
         health_mean.append(healthy_val[int(spli[1])-1])  ## collecting the healthy equivalant 6MWD
         
         slope_fitbit.append(m1)     ## collectig the slope of Fitbit 6MWD
         
         
   




