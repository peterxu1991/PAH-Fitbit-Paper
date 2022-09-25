# -*- coding: utf-8 -*-
"""
Created on Tue Jul  5 10:39:20 2022

@author: Peter Xu

"""
import pandas as pd
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt

## read all the data in each sheet from excel file
df_Fitbit = pd.read_excel(r'C:\Users\zxu11\OneDrive\Desktop\JHU Paper\sensi_analysis.xlsx',sheet_name = 0)
df_clinic = pd.read_excel(r'C:\Users\zxu11\OneDrive\Desktop\JHU Paper\sensi_analysis.xlsx',sheet_name = 1)

### Performer vs Non-performer

df1 = df_clinic.iloc[:,1:] #all clinical variables

p_pm = []

for col_no,col in enumerate(df1):
    
    grp_high = df1.iloc[np.r_[0,2,4,8,9,12:17,18,19],col_no]
    grp_low = df1.iloc[np.r_[1,3,5:8,10,11,17,20,21],col_no]
    
    tStat, pValue = stats.mannwhitneyu(grp_low,grp_high,nan_policy='omit') #run independent sample T-Test
    p_pm.append(pValue)
    
    if pValue < 0.05:
        clin_raw = pd.concat([grp_low,grp_high],axis=1,ignore_index=True)
        clin_raw.columns = ['Under_Performer','Performer']
        with pd.ExcelWriter(r'C:\Users\zxu11\OneDrive\Desktop\JHU Paper\clinic_group\Performance\%s.xlsx'%(col)) as writer:
            clin_raw.to_excel(writer)
   
# p_pm = pd.DataFrame(p_pm)
# p_pm.index = df1.columns

# p_pm.to_excel(r'C:\Users\zxu11\OneDrive\Desktop\JHU Paper\performance.xlsx')

# ### RHR vs all clinic

# df_RHR = df_Fitbit['HR(SC=0)mean'] # average daily steps
# df1 = df_clinic.iloc[:,1:] #all clinical variables

# RHR_range = range(70,92,2)
# p_RHR = []

# for col_no,col in enumerate(df1):
    
#     df_grp = pd.concat([df_RHR,df1[col]],axis=1,ignore_index=True) 
#     df_grp.columns = ['RHR',col]
#     p_RHR.append([])
    
#     for RHR in RHR_range:
#         grp_low = df_grp[col][df_grp['RHR']<=RHR]
#         grp_high = df_grp[col][df_grp['RHR']>RHR]
        
#         if (len(grp_low) >= 5) & (len(grp_high) >= 5):
#             tStat, pValue = stats.mannwhitneyu(grp_low,grp_high,nan_policy='omit') #run independent sample T-Test
#             p_RHR[col_no].append(pValue)
            
#             if (RHR == 82) & (pValue < 0.05):
#                 clin_raw = pd.concat([grp_low,grp_high],axis=1,ignore_index=True)
#                 clin_raw.columns = ['Fitbit_lower','Fitbit_higher']
#                 with pd.ExcelWriter(r'C:\Users\zxu11\OneDrive\Desktop\JHU Paper\clinic_group\RHR\%s.xlsx'%(col)) as writer:
#                     clin_raw.to_excel(writer)
#         else:
#             pValue = float('NaN')
#             p_RHR[col_no].append(pValue)
        
            
# p_RHR = pd.DataFrame(p_RHR)
# p_RHR.index = df1.columns
# p_RHR = p_RHR.T
# list_RHR = list(RHR_range)
# p_RHR.index = [str(x) for x in list_RHR]

# ## plot the number of significant clinical vars versus threshold values
# temp_p = p_RHR[p_RHR < 0.05]
# cnt_RHR = temp_p.count(axis = 1)

# df_cut = pd.DataFrame(list_RHR)
# cnt_RHR.index = list(range(0,len(list_RHR)))
# sigP_thres = pd.concat([df_cut,cnt_RHR],axis=1)
# sigP_thres.columns = ['RHR_thresholds','Number of Significant pValue']
# sigP_thres.to_excel(r'C:\Users\zxu11\OneDrive\Desktop\JHU Paper\SigP vs Cutoff\RHR.xlsx')

# fig1 = plt.figure(figsize = (15,15))
# plt.plot(list_RHR, cnt_RHR,lw=2) 
# plt.xlabel(df_Fitbit.columns[1])
# plt.ylabel('Significant Clinical Param No.')
# plt.title('No of Significant Parameters vs %s Cut-offs'%(df_Fitbit.columns[1]))
# fig1.savefig(r'C:\Users\zxu11\OneDrive\Desktop\JHU Paper\plots\RHR.jpg')



# ### RHR_sk vs all clinic

# df_HRsk = df_Fitbit['HR(SC=0)sk'] # average daily steps
# df1 = df_clinic.iloc[:,1:] #all clinical variables

# HRsk_range = np.arange(0.2,1.6,0.1)
# p_HRsk = []

# for col_no,col in enumerate(df1):
    
#     df_grp = pd.concat([df_HRsk,df1[col]],axis=1,ignore_index=True) 
#     df_grp.columns = ['HRsk',col]
#     p_HRsk.append([])
    
#     for HRsk in HRsk_range:
#         grp_low = df_grp[col][df_grp['HRsk']<=HRsk]
#         grp_high = df_grp[col][df_grp['HRsk']>HRsk]
        
#         if (len(grp_low) >= 5) & (len(grp_high) >= 5):
#             tStat, pValue = stats.mannwhitneyu(grp_low,grp_high,nan_policy='omit') #run independent sample T-Test
#             p_HRsk[col_no].append(pValue)
            
#             if ((HRsk > 0.9) & (HRsk <1.1)) & (pValue < 0.05):
#                 clin_raw1 = pd.concat([grp_low,grp_high],axis=1,ignore_index=True)
#                 clin_raw1.columns = ['Fitbit_lower','Fitbit_higher']
#                 with pd.ExcelWriter(r'C:\Users\zxu11\OneDrive\Desktop\JHU Paper\clinic_group\HRsk\%s.xlsx'%(col)) as writer:
#                     clin_raw1.to_excel(writer)
#         else:
#             pValue = float('NaN')
#             p_HRsk[col_no].append(pValue)
        
# p_HRsk = pd.DataFrame(p_HRsk)
# p_HRsk.index = df1.columns
# p_HRsk=p_HRsk.T
# list_HRsk = list(HRsk_range)
# p_HRsk.index = [str(x) for x in list_HRsk]

# ## plot the number of significant clinical vars versus threshold values
# temp_p = p_HRsk[p_HRsk < 0.05]
# cnt_HRsk = temp_p.count(axis = 1)

# df_cut = pd.DataFrame(list_HRsk)
# cnt_HRsk.index = list(range(0,len(list_HRsk)))
# sigP_thres = pd.concat([df_cut,cnt_HRsk],axis=1)
# sigP_thres.columns = ['HRsk_thresholds','Number of Significant pValue']
# sigP_thres.to_excel(r'C:\Users\zxu11\OneDrive\Desktop\JHU Paper\SigP vs Cutoff\HRsk.xlsx')

# fig1 = plt.figure(figsize = (15,15))
# plt.plot(list_HRsk, cnt_HRsk,lw=2) 
# plt.xlabel(df_Fitbit.columns[2])
# plt.ylabel('Significant Clinical Param No.')
# plt.title('No of Significant Parameters vs %s Cut-offs'%(df_Fitbit.columns[2]))
# fig1.savefig(r'C:\Users\zxu11\OneDrive\Desktop\JHU Paper\plots\HRsk.jpg')

# ### HR_SC>0 vs all clinic

# df_SC0 = df_Fitbit['HR(SC>0)mean'] # average daily steps
# df1 = df_clinic.iloc[:,1:] #all clinical variables

# SC0_range = range(80,110,5)
# p_SC0 = []

# for col_no,col in enumerate(df1):
    
#     df_grp = pd.concat([df_SC0,df1[col]],axis=1,ignore_index=True) 
#     df_grp.columns = ['SC0',col]
#     p_SC0.append([])
    
#     for SC0 in SC0_range:
#         grp_low = df_grp[col][df_grp['SC0']<=SC0]
#         grp_high = df_grp[col][df_grp['SC0']>SC0]
        
#         if (len(grp_low) >= 5) & (len(grp_high) >= 5):
#             tStat, pValue = stats.mannwhitneyu(grp_low,grp_high,nan_policy='omit') #run independent sample T-Test
#             p_SC0[col_no].append(pValue)
            
#             if (SC0 == 95) & (pValue < 0.05):
#                 clin_raw2 = pd.concat([grp_low,grp_high],axis=1,ignore_index=True)
#                 clin_raw2.columns = ['Fitbit_lower','Fitbit_higher']
#                 with pd.ExcelWriter(r'C:\Users\zxu11\OneDrive\Desktop\JHU Paper\clinic_group\HRSC1\%s.xlsx'%(col)) as writer:
#                     clin_raw2.to_excel(writer)
#         else:
#             pValue = float('NaN')
#             p_SC0[col_no].append(pValue)
        
# p_SC0 = pd.DataFrame(p_SC0)
# p_SC0.index = df1.columns
# p_SC0=p_SC0.T
# list_SC0 = list(SC0_range)
# p_SC0.index = [str(x) for x in list_SC0]

# ## plot the number of significant clinical vars versus threshold values
# temp_p = p_SC0[p_SC0 < 0.05]
# cnt_SC0 = temp_p.count(axis = 1)

# df_cut = pd.DataFrame(list_SC0)
# cnt_SC0.index = list(range(0,len(list_SC0)))
# sigP_thres = pd.concat([df_cut,cnt_SC0],axis=1)
# sigP_thres.columns = ['HRSC1_thresholds','Number of Significant pValue']
# sigP_thres.to_excel(r'C:\Users\zxu11\OneDrive\Desktop\JHU Paper\SigP vs Cutoff\HRSC1.xlsx')

# fig1 = plt.figure(figsize = (15,15))
# plt.plot(list_SC0, cnt_SC0,lw=2) 
# plt.xlabel(df_Fitbit.columns[3])
# plt.ylabel('Significant Clinical Param No.')
# plt.title('No of Significant Parameters vs %s Cut-offs'%(df_Fitbit.columns[3]))
# fig1.savefig(r'C:\Users\zxu11\OneDrive\Desktop\JHU Paper\plots\HR_SC0.jpg')

# ### average_daily_steps vs all clinic

# df_ave_steps = df_Fitbit['avg_daily_steps'] # average daily steps
# df1 = df_clinic.iloc[:,1:] #all clinical variables

# step_range = range(1500,10000,500)
# p_steps = []

# for col_no,col in enumerate(df1):
    
#     df_grp = pd.concat([df_ave_steps,df1[col]],axis=1,ignore_index=True) 
#     df_grp.columns = ['avg_steps',col]
#     p_steps.append([])
    
#     for step in step_range:
#         grp_low = df_grp[col][df_grp['avg_steps']<=step]
#         grp_high = df_grp[col][df_grp['avg_steps']>step]
        
#         if (len(grp_low) >= 5) & (len(grp_high) >= 5):
#             tStat, pValue = stats.mannwhitneyu(grp_low,grp_high,nan_policy='omit') #run independent sample T-Test
#             p_steps[col_no].append(pValue)
            
#             if (step == 5000) & (pValue < 0.05):
#                 clin_raw3 = pd.concat([grp_low,grp_high],axis=1,ignore_index=True)
#                 clin_raw3.columns = ['Fitbit_lower','Fitbit_higher']
#                 with pd.ExcelWriter(r'C:\Users\zxu11\OneDrive\Desktop\JHU Paper\clinic_group\daily_steps\%s.xlsx'%(col)) as writer:
#                     clin_raw3.to_excel(writer)
#         else:
#             pValue = float('NaN')
#             p_steps[col_no].append(pValue)
        
# p_steps = pd.DataFrame(p_steps)
# p_steps.index = df1.columns
# p_steps=p_steps.T
# list_step = list(step_range)
# p_steps.index = [str(x) for x in list_step]

# ## plot the number of significant clinical vars versus threshold values
# temp_p = p_steps[p_steps < 0.05]
# cnt_step = temp_p.count(axis = 1)

# df_cut = pd.DataFrame(list_step)
# cnt_step.index = list(range(0,len(list_step)))
# sigP_thres = pd.concat([df_cut,cnt_step],axis=1)
# sigP_thres.columns = ['Daily_steps_thresholds','Number of Significant pValue']
# sigP_thres.to_excel(r'C:\Users\zxu11\OneDrive\Desktop\JHU Paper\SigP vs Cutoff\daily_steps.xlsx')

# fig1 = plt.figure(figsize = (15,15))
# plt.plot(list_step, cnt_step,lw=2) 
# plt.xlabel(df_Fitbit.columns[4])
# plt.ylabel('Significant Clinical Param No.')
# plt.title('No of Significant Parameters vs %s Cut-offs'%(df_Fitbit.columns[4]))
# fig1.savefig(r'C:\Users\zxu11\OneDrive\Desktop\JHU Paper\plots\%s.jpg'%(df_Fitbit.columns[4]))

# ### fitness_slope vs all clinic

# df_fit_slope = df_Fitbit['fitness_slope'] # average daily steps
# df1 = df_clinic.iloc[:,1:] #all clinical variables

# fit_slope_range = np.arange(0.03,0.3,0.03)
# p_fit_slope = []

# for col_no,col in enumerate(df1):
    
#     df_grp = pd.concat([df_fit_slope,df1[col]],axis=1,ignore_index=True) 
#     df_grp.columns = ['fit_slope',col]
#     p_fit_slope.append([])
    
#     for fit_slope in fit_slope_range:
#         grp_low = df_grp[col][df_grp['fit_slope'] <= fit_slope]
#         grp_high = df_grp[col][df_grp['fit_slope'] > fit_slope]
        
#         if (len(grp_low) >= 5) & (len(grp_high) >= 5):
#             tStat, pValue = stats.mannwhitneyu(grp_low,grp_high,nan_policy='omit') #run independent sample T-Test
#             p_fit_slope[col_no].append(pValue)
            
#             if (fit_slope == 0.15) & (pValue < 0.05):
#                 clin_raw4 = pd.concat([grp_low,grp_high],axis=1,ignore_index=True)
#                 clin_raw4.columns = ['Fitbit_lower','Fitbit_higher']
#                 with pd.ExcelWriter(r'C:\Users\zxu11\OneDrive\Desktop\JHU Paper\clinic_group\fitness\%s.xlsx'%(col)) as writer:
#                     clin_raw4.to_excel(writer)
#         else:
#             pValue = float('NaN')
#             p_fit_slope[col_no].append(pValue)
        
# p_fit_slope = pd.DataFrame(p_fit_slope)
# p_fit_slope.index = df1.columns
# p_fit_slope=p_fit_slope.T
# list_fit_slope = list(fit_slope_range)
# p_fit_slope.index = [str(x) for x in list_fit_slope]

# ## plot the number of significant clinical vars versus threshold values
# temp_p = p_fit_slope[p_fit_slope < 0.05]
# cnt_fit_slope = temp_p.count(axis = 1)

# df_cut = pd.DataFrame(list_fit_slope)
# cnt_fit_slope.index = list(range(0,len(list_fit_slope)))
# sigP_thres = pd.concat([df_cut,cnt_fit_slope],axis=1)
# sigP_thres.columns = ['fitness_thresholds','Number of Significant pValue']
# sigP_thres.to_excel(r'C:\Users\zxu11\OneDrive\Desktop\JHU Paper\SigP vs Cutoff\fitness.xlsx')

# fig1 = plt.figure(figsize = (15,15))
# plt.plot(list_fit_slope, cnt_fit_slope,lw=2) 
# plt.xlabel(df_Fitbit.columns[5])
# plt.ylabel('Significant Clinical Param No.')
# plt.title('No of Significant Parameters vs %s Cut-offs'%(df_Fitbit.columns[5]))
# fig1.savefig(r'C:\Users\zxu11\OneDrive\Desktop\JHU Paper\plots\%s.jpg'%(df_Fitbit.columns[5]))

# ### ambulation_product vs all clinic

# df_ambu = df_Fitbit['ambulation_product'] # average daily steps
# df1 = df_clinic.iloc[:,1:] #all clinical variables

# ambu_range = range(500,1500,100)
# p_ambu = []

# for col_no,col in enumerate(df1):
    
#     df_grp = pd.concat([df_ambu,df1[col]],axis=1,ignore_index=True) 
#     df_grp.columns = ['ambu',col]
#     p_ambu.append([])
    
#     for ambu in ambu_range:
#         grp_low = df_grp[col][df_grp['ambu'] <= ambu]
#         grp_high = df_grp[col][df_grp['ambu'] > ambu]
        
#         if (len(grp_low) >= 5) & (len(grp_high) >= 5):
#             tStat, pValue = stats.mannwhitneyu(grp_low,grp_high,nan_policy='omit') #run independent sample T-Test
#             p_ambu[col_no].append(pValue)
            
#             if (ambu == 1000) & (pValue < 0.05):
#                 clin_raw5 = pd.concat([grp_low,grp_high],axis=1,ignore_index=True)
#                 clin_raw5.columns = ['Fitbit_lower','Fitbit_higher']
#                 with pd.ExcelWriter(r'C:\Users\zxu11\OneDrive\Desktop\JHU Paper\clinic_group\Ambulation\%s.xlsx'%(col)) as writer:
#                     clin_raw5.to_excel(writer)
#         else:
#             pValue = float('NaN')
#             p_ambu[col_no].append(pValue)
        
# p_ambu = pd.DataFrame(p_ambu)
# p_ambu.index = df1.columns
# p_ambu = p_ambu.T
# list_ambu = list(ambu_range)
# p_ambu.index = [str(x) for x in list_ambu]

# ## plot the number of significant clinical vars versus threshold values
# temp_p = p_ambu[p_ambu < 0.05]
# cnt_ambu = temp_p.count(axis = 1)

# df_cut = pd.DataFrame(list_ambu)
# cnt_ambu.index = list(range(0,len(list_ambu)))
# sigP_thres = pd.concat([df_cut,cnt_ambu],axis=1)
# sigP_thres.columns = ['Ambulation_thresholds','Number of Significant pValue']
# sigP_thres.to_excel(r'C:\Users\zxu11\OneDrive\Desktop\JHU Paper\SigP vs Cutoff\Ambulation.xlsx')

# fig1 = plt.figure(figsize = (15,15))
# plt.plot(list_ambu, cnt_ambu,lw=2) 
# plt.xlabel(df_Fitbit.columns[6])
# plt.ylabel('Significant Clinical Param No.')
# plt.title('No of Significant Parameters vs %s Cut-offs'%(df_Fitbit.columns[6]))
# fig1.savefig(r'C:\Users\zxu11\OneDrive\Desktop\JHU Paper\plots\%s.jpg'%(df_Fitbit.columns[6]))

# ### usage vs all clinic

# df_usage = df_Fitbit['usage'] # average daily steps
# df1 = df_clinic.iloc[:,1:] #all clinical variables

# usage_range = np.arange(0.8,0.98,0.02)
# p_usage = []

# for col_no,col in enumerate(df1):
    
#     df_grp = pd.concat([df_usage,df1[col]],axis=1,ignore_index=True) 
#     df_grp.columns = ['usage',col]
#     p_usage.append([])
    
#     for usage in usage_range:
#         grp_low = df_grp[col][df_grp['usage']<=usage]
#         grp_high = df_grp[col][df_grp['usage']>usage]
        
#         if (len(grp_low) >= 5) & (len(grp_high) >= 5):
#             tStat, pValue = stats.mannwhitneyu(grp_low,grp_high,nan_policy='omit') #run independent sample T-Test
#             p_usage[col_no].append(pValue)
            
#             if ((usage > 0.93) & (usage < 0.95)) & (pValue < 0.05):
#                 clin_raw6 = pd.concat([grp_low,grp_high],axis=1,ignore_index=True)
#                 clin_raw6.columns = ['Fitbit_lower','Fitbit_higher']
#                 with pd.ExcelWriter(r'C:\Users\zxu11\OneDrive\Desktop\JHU Paper\clinic_group\usage\%s.xlsx'%(col)) as writer:
#                     clin_raw6.to_excel(writer)
#         else:
#             pValue = float('NaN')
#             p_usage[col_no].append(pValue)
        
# p_usage = pd.DataFrame(p_usage)
# p_usage.index = df1.columns
# p_usage=p_usage.T
# list_usage = list(usage_range)
# p_usage.index = [str(x) for x in list_usage]

# ## plot the number of significant clinical vars versus threshold values
# temp_p = p_usage[p_usage < 0.05]
# cnt_usage = temp_p.count(axis = 1)

# df_cut = pd.DataFrame(list_usage)
# cnt_usage.index = list(range(0,len(list_usage)))
# sigP_thres = pd.concat([df_cut,cnt_usage],axis=1)
# sigP_thres.columns = ['usage_thresholds','Number of Significant pValue']
# sigP_thres.to_excel(r'C:\Users\zxu11\OneDrive\Desktop\JHU Paper\SigP vs Cutoff\usage.xlsx')

# fig1 = plt.figure(figsize = (15,15))
# plt.plot(list_usage, cnt_usage,lw=2) 
# plt.xlabel(df_Fitbit.columns[7])
# plt.ylabel('Significant Clinical Param No.')
# plt.title('No of Significant Parameters vs %s Cut-offs'%(df_Fitbit.columns[7]))
# fig1.savefig(r'C:\Users\zxu11\OneDrive\Desktop\JHU Paper\plots\%s.jpg'%(df_Fitbit.columns[7]))

### FL6MWT vs all clinic

df_MWT = df_Fitbit['FL6MWT'] # average daily steps
df1 = df_clinic.iloc[:,1:] #all clinical variables

MWT_range = range(200,450,20)
p_MWT = []

for col_no,col in enumerate(df1):
    
    df_grp = pd.concat([df_MWT,df1[col]],axis=1,ignore_index=True) 
    df_grp.columns = ['MWT',col]
    p_MWT.append([])
    
    for MWT in MWT_range:
        grp_low = df_grp[col][df_grp['MWT']<=MWT]
        grp_high = df_grp[col][df_grp['MWT']>MWT]
        
        if (len(grp_low) >= 5) & (len(grp_high) >= 5):
            tStat, pValue = stats.mannwhitneyu(grp_low,grp_high,nan_policy='omit') #run independent sample T-Test
            p_MWT[col_no].append(pValue)
            
            if (MWT == 400) & (pValue < 0.05):
                clin_raw7 = pd.concat([grp_low,grp_high],axis=1,ignore_index=True)
                clin_raw7.columns = ['Fitbit_lower','Fitbit_higher']
                with pd.ExcelWriter(r'C:\Users\zxu11\OneDrive\Desktop\JHU Paper\clinic_group\6MWT_400m\%s.xlsx'%(col)) as writer:
                    clin_raw7.to_excel(writer)
        else:
            pValue = float('NaN')
            p_MWT[col_no].append(pValue)
        
p_MWT = pd.DataFrame(p_MWT)
p_MWT.index = df1.columns
p_MWT = p_MWT.T
list_MWT = list(MWT_range)
p_MWT.index = [str(x) for x in list_MWT]

## plot the number of significant clinical vars versus threshold values
temp_p = p_MWT[p_MWT < 0.05]
cnt_MWT = temp_p.count(axis = 1)

df_cut = pd.DataFrame(list_MWT)
cnt_MWT.index = list(range(0,len(list_MWT)))
sigP_thres = pd.concat([df_cut,cnt_MWT],axis=1)
sigP_thres.columns = ['6MWT_thresholds','Number of Significant pValue']
sigP_thres.to_excel(r'C:\Users\zxu11\OneDrive\Desktop\JHU Paper\SigP vs Cutoff\6MWT.xlsx')

fig1 = plt.figure(figsize = (15,15))
plt.plot(list_MWT, cnt_MWT,lw=2) 
plt.xlabel(df_Fitbit.columns[8])
plt.ylabel('Significant Clinical Param No.')
plt.title('No of Significant Parameters vs %s Cut-offs'%(df_Fitbit.columns[8]))
fig1.savefig(r'C:\Users\zxu11\OneDrive\Desktop\JHU Paper\plots\%s.jpg'%(df_Fitbit.columns[8]))

### 6MWT_slope vs all clinic

# =============================================================================
# df_MWslope = df_Fitbit['6MW_slope'] # average daily steps
# df1 = df_clinic.iloc[:,1:] #all clinical variables
# 
# MWslope_range = range(-10,5,1)
# p_MWslope = []
# 
# for col_no,col in enumerate(df1):
#     
#     df_grp = pd.concat([df_MWslope,df1[col]],axis=1,ignore_index=True) 
#     df_grp.columns = ['MWslope',col]
#     p_MWslope.append([])
#     
#     for MWslope in MWslope_range:
#         grp_low = df_grp[col][df_grp['MWslope']<=MWslope]
#         grp_high = df_grp[col][df_grp['MWslope']>MWslope]
#         
#         if (len(grp_low) > 5) & (len(grp_high) > 5):
#             tStat, pValue = stats.mannwhitneyu(grp_low,grp_high) #run independent sample T-Test
#             p_MWslope[col_no].append(pValue)
#         else:
#             pValue = float('NaN')
#             p_MWslope[col_no].append(pValue)
#         
# p_MWslope = pd.DataFrame(p_MWslope)
# p_MWslope.index = df1.columns
# p_MWslope=p_MWslope.T
# list_MWslope = list(MWslope_range)
# p_MWslope.index = [str(x) for x in list_MWslope]
# 
# ## plot the number of significant clinical vars versus threshold values
# temp_p = p_MWslope[p_MWslope < 0.05]
# cnt_MWslope = temp_p.count(axis = 1)
# 
# fig1 = plt.figure(figsize = (15,15))
# plt.plot(list_MWslope, cnt_MWslope,lw=2) 
# plt.xlabel(df_Fitbit.columns[9])
# plt.ylabel('Significant Clinical Param No.')
# plt.title('No of Significant Parameters vs %s Cut-offs'%(df_Fitbit.columns[9]))
# fig1.savefig(r'C:\Users\zxu11\OneDrive\Desktop\JHU Paper\plots\%s.jpg'%(df_Fitbit.columns[9]))
# =============================================================================

### save all the data in one excel

# # create a excel writer object

# with pd.ExcelWriter(r'C:\Users\zxu11\OneDrive\Desktop\JHU Paper\sens_results_mannWhitney_v2.xlsx') as writer:
#     p_RHR.to_excel(writer,sheet_name = 'RHR')
#     p_HRsk.to_excel(writer,sheet_name = 'HR(SC=0)sk')
#     p_SC0.to_excel(writer,sheet_name = 'HR(SC>0)')
#     p_steps.to_excel(writer,sheet_name = 'Avg_daily_steps')
#     p_fit_slope.to_excel(writer,sheet_name = 'fitness_slope')
#     p_ambu.to_excel(writer,sheet_name = 'Ambulation')
#     p_usage.to_excel(writer,sheet_name = 'Usage')
#     p_MWT.to_excel(writer,sheet_name = 'FL6MWT')



