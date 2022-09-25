import numpy as np
import pandas as pd
import math
import numpy.matlib
from sklearn.cluster import KMeans

# =============================================================================
#  Make sure to have 'PAH_6MWD2.xlsx' file in the same directory
# =============================================================================
xls = pd.ExcelFile('PAH_6MWD2.xlsx')
sheet_names = xls.sheet_names

MM = []
SD = []
Slope = []
Pred = []
pred_2 =[]
predict =[]
Actual = []
TT =0
MMM =[]
for i in range(len(sheet_names)): 


        # =============================================================================
        #   Reading all the Exel File Sheets
        # =============================================================================
        xls = pd.ExcelFile('PAH_6MWD2.xlsx')
        sheet_names = xls.sheet_names

        df = pd.read_excel(xls, sheet_names[i])
        
        First_val = df.iloc[0,1]
        Second_val = df.iloc[-1,1]
        
        if  not (math.isnan(Second_val) or math.isnan(First_val)):
            
            FitBit_arr = np.array(df.iloc[0+1: -1,1+1]).reshape(-1,1)
            arr=FitBit_arr[np.logical_not(np.isnan(FitBit_arr))]
            
            
# =============================================================================
#             Calculating the disatcen to the clinic line
# =============================================================================
            TT= np.vstack((TT,arr.reshape(-1,1)))
            m0 = (Second_val-First_val)/(df.shape[0]-1)
            m, b = np.polyfit([1,df.shape[0]],[First_val, Second_val], 1)
            m1, b1 = np.polyfit(range(2,len(arr)+2),arr[:], 1)
                        
            pp = Second_val - (m1*df.shape[0]+b1)            
            Pred.append(pp)
            pred_2.append((m1*df.shape[0]+b1))
            
            
            predict.append((m1*df.shape[0]+b1))
            Actual.append(Second_val)
            
            Dist = []
            
            MMM.append(np.mean(arr))
            for j in range(2,len(arr)+2): 
                p3 = np.array([j,arr[j-2]])   
                d = p3[1]-(m*p3[0]+b)
                Dist.append(d)
             
            Dist = np.vstack(Dist)
            
            MU = np.mean(Dist)
            MM.append(MU)
      
            sd = np.std(Dist)
            SD.append(sd)
Data = np.hstack((np.vstack(MM),np.vstack(SD)))        


X = Data[:,0].reshape(-1,1)
X = np.round(X,2)
kmeans = KMeans(n_clusters=2, random_state = 0).fit(X)   





sheet_names0 =sheet_names
sheet_names.remove('PAH8')     ## PAH 8 is not included for the analysis 


A = kmeans.labels_.reshape(-1,1)


kmeans.labels_ = (-(kmeans.labels_ -1))  
kmeans.labels_[kmeans.labels_ == 0] =2
Data  = np.hstack((Data,kmeans.labels_.reshape(-1,1)))        
Data0 =  np.vstack(Pred)

New = Data0[kmeans.labels_ == 1]  
df_T =  pd.DataFrame(Data,columns =['Mean','SD','Group'])  ## Final table with value input to the k-means and Labels

## Group 1 : Performers   &&  Group 2: Underperformers