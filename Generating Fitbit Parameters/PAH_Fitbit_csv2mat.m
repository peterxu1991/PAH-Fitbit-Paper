%%% convert all the csv file to the mat file
clear;clc;close all;
%% read data
PathName = 'C:\Users\Peter Xu\OneDrive - Johns Hopkins University\JHU\Projects\Fitbit_IRB_28FEB2018\PAH_Fitbit\Minute_Data_all\';
date_range = '_20160101_20170822';

for i = 1 : 30
    
    empty_sub = [6,16,24];
    
    if ~isempty(find(empty_sub == i))  
        continue;
    else
        subi_name = ['FBPH' num2str(i,'%02d')];
        File_format = '.csv';
        
        File_HR = [subi_name '_heartrate_1min' date_range File_format];
        File_SC = [subi_name '_minuteStepsNarrow' date_range File_format];
        File_ACT = [subi_name '_minuteIntensitiesNarrow' date_range File_format];
        % heart rate and its time stamp
        RawHR = ReadPAHFitbit(PathName,File_HR);
        % step count and its time stamp
        RawSC = ReadPAHFitbit(PathName,File_SC);
        % Activity level and its time stamp
        RawACT = ReadPAHFitbit(PathName,File_ACT);
        
        temp1 = RawHR.value;
        temp2 = RawSC.value;
        temp3 = RawACT.value;
        time_HR = RawHR.time;
        time_SC = RawSC.time;
        
        [subi_time,IA,IB] = intersect(time_HR,time_SC);
        
        subi_hr = temp1(IA);
        subi_step = temp2(IB);
        
        subi_act = temp3(IA);
        
        RawFBData(i).hr = subi_hr;
        RawFBData(i).step = subi_step;
        RawFBData(i).act = subi_act;
        RawFBData(i).time = subi_time;
        RawFBData(i).name = subi_name;
        
        
    end
end
% % save data to mat file
% save 'PAH_Fitbit.mat' RawFBData;