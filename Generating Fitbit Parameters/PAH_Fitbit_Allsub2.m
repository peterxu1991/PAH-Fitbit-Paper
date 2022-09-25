%%% Fitbit step & Heart rate

clear;clc;close all;
% read the age list of each subject
PathName_age = 'C:\Users\zxu11\OneDrive - Johns Hopkins University\JHU\Projects\Fitbit_IRB_28FEB2018\PAH_Fitbit\';
FileName_age = 'PAH_Fitbit_age.xlsx';
[temp_num,~,~] = xlsread([PathName_age,FileName_age]);
sub_ID = temp_num(:,1);
[~,sort_ind] = sort(sub_ID);
age_list = round(temp_num(sort_ind,3));

FileName_date = 'dates.xlsx';
[temp_num2,txt2,~] = xlsread([PathName_age,FileName_date]);
txt2{8,2} = NaN;
txt2{9,3} = NaN;txt2{19,3} = NaN;txt2{20,3} = NaN;txt2{30,3} = NaN;txt2{31,3} = NaN;

start_date = zeros(30,1);
end_date = zeros(30,1);

for i = 1 : length(temp_num2)
    
    start_date(i) = datenum(txt2{i+1,2});
    end_date(i) = datenum(txt2{i+1,3});
    
end

%% read data

PathName = 'C:\Users\zxu11\OneDrive - Johns Hopkins University\JHU\Projects\Fitbit_IRB_28FEB2018\PAH_Fitbit\Minute_Data_mat\';
FileName = 'PAH_Fitbit.mat';
load([PathName,FileName]);

SleepFile = 'PAH_SleepData.mat';
load([PathName,SleepFile]);
%% pre-allocation of wear pct daily

base_week = zeros(length(RawFBData),3);
step_avg_hour = zeros(24,2,length(RawFBData));
hr_avg_hour = zeros(24,2,length(RawFBData));
act_avg_hour = zeros(24,2,length(RawFBData));

for n = 1 : length(RawFBData)
    
    empty_sub = [6,16,18,24,25,29];
    
    if ~isempty(find(empty_sub == n))
        continue;
    else
        
        subi_step = RawFBData(n).step;
        subi_hr = RawFBData(n).hr;
        subi_time = RawFBData(n).time;
        subi_name = RawFBData(n).name;
        subi_act = RawFBData(n).act;
        subi_age = age_list(n);
        subi_sleep = RawFBSleep(n).sleep;
        
        % remove the repeated rows with same time
        [subi_time,ia, ~] = unique(subi_time);
        subi_hr = subi_hr(ia);
        subi_step = subi_step(ia);
        subi_act = subi_act(ia);
        
        % age-predicted maximum heart rate (refer to Tanaka et al. 2001)
        ap_mhr = 208 - 0.7 * subi_age;
        
        L = 1: 1 : length(subi_time);
        
        % truncate the data based on the start and end date
        if ~isnan(start_date(n)) && ~isnan(end_date(n))
            ind_start = find(subi_time > start_date(n));
            ind_end = find(subi_time < end_date(n));
            
            ind = intersect(ind_start,ind_end);
            
            subi_hr = subi_hr(ind);
            subi_step = subi_step(ind);
            subi_act = subi_act(ind);
            subi_time = subi_time(ind);
        end
        
        % convert time with datenum format to datetime format
        t = datetime(unique(subi_time),'ConvertFrom','datenum');
        
        %% plot the wear time fraction for each day
        
        %         plot_time_fract(t,subi_act,subi_name);
        
        %% Process the sleep data
        %% pre-process data
        % remove zero heart rate points
        
        [rx,~] = find(subi_hr ==0);
        rm_hr0 = length(rx);
        subi_hr(rx) = [];
        subi_step(rx) = [];
        subi_act(rx)=[];
        t(rx) = [];
        
        % keep heart rate data if it is between 20 and age-predicted max heart rate + SD (10bpm) (refer to Avram R et al._2019)
        [rx1,~] = find(subi_hr > ap_mhr + 10 | subi_hr < 20);
        rm_hrout = length(rx1);
        subi_hr(rx1) = [];
        subi_step(rx1) = [];
        subi_act(rx1) = [];
        t(rx1) = [];
        
        %     t.Format = 'dd-MMM-yyyy';
        %% plot the missing data distribution
        [week_no,compliance, missing_hours,day_pct_week,night_pct_week,area_step_week,max_hour_step, ...
            core_hour1, core_hour2, step_en, hr_en, r_step_en, r_hr_en] = Fitbit_MissingData2(subi_hr,subi_step,subi_act,t,subi_name,ap_mhr);
        %
        %% process the data on the daily basis at different window size in the temporal space
        [act_abs,act_pct,fa, hr_mu0,hr_sig0,hr_sk0,hr_ku0,hr_ks0,hr_mu1,hr_sig1,...
            hr_sk1,hr_ku1,hr_ks1,sc_mu,sc_sig,sc_sk,sc_ku,sc_ks,slope,intercept,area1,...
            max_hr_week,hr_max_sc,max_step_week,sc_max_hr,avg_top_sc,avg_top_hr,...
            avg_6MWT_sc,avg_6MWT_hr,Pred_SC170,Pred_HR170,weekly_raw] = Fitbit_DataProcess4(subi_hr,subi_step,subi_act,t,subi_name);
        
        %% save the variables to mat file
        daily_PAH(n).age = subi_age;
        daily_PAH(n).name = subi_name;
        
        daily_PAH(n).hr_mu0 = hr_mu0;
        daily_PAH(n).hr_sig0 = hr_sig0;
        daily_PAH(n).hr_sk0 = hr_sk0;
        daily_PAH(n).hr_ku0 = hr_ku0;
        daily_PAH(n).hr_ks0 = hr_ks0;
        
        daily_PAH(n).hr_mu1 = hr_mu1;
        daily_PAH(n).hr_sig1 = hr_sig1;
        daily_PAH(n).hr_sk1 = hr_sk1;
        daily_PAH(n).hr_ku1 = hr_ku1;
        daily_PAH(n).hr_ks1 = hr_ks1;
        
        daily_PAH(n).sc_mu = sc_mu;
        daily_PAH(n).sc_sig = sc_sig;
        daily_PAH(n).sc_sk = sc_sk;
        daily_PAH(n).sc_ku = sc_ku;
        daily_PAH(n).sc_ks = sc_ks;
        
        daily_PAH(n).slope = slope;
        daily_PAH(n).intercept = intercept;
        daily_PAH(n).area = area1;
        
        daily_PAH(n).act_abs = act_abs;
        daily_PAH(n).act_pct = act_pct;
        
        
        daily_PAH(n).facti = fa;
        daily_PAH(n).day_pct = day_pct_week;
        daily_PAH(n).night_pct = night_pct_week;
        daily_PAH(n).compliance = compliance;
        daily_PAH(n).missing_hours = missing_hours;
        
        daily_PAH(n).area_step_week = area_step_week;
        daily_PAH(n).max_hour_step = max_hour_step;
        
        daily_PAH(n).core_hour1 = core_hour1;
        
        daily_PAH(n).max_step = max_step_week;
        daily_PAH(n).hr_max_sc = hr_max_sc;
        daily_PAH(n).max_hr = max_hr_week;
        daily_PAH(n).sc_max_hr = sc_max_hr;
        
        daily_PAH(n).mean_top_sc = avg_top_sc;
        daily_PAH(n).mean_top_hr = avg_top_hr;
        daily_PAH(n).MWT_sc = avg_6MWT_sc;
        daily_PAH(n).MWT_hr = avg_6MWT_hr;
        
        daily_PAH(n).SC170 = Pred_SC170;
        daily_PAH(n).Hr170 = Pred_HR170;
        daily_PAH(n).weekly_raw = weekly_raw;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %                 daily_PAH(n).core_hour2 = core_hour2;
        %
        %                 daily_PAH(n).step_en = step_en;
        %                 daily_PAH(n).hr_en = hr_en;
        %
        %                 daily_PAH(n).r_step_en = r_step_en;
        %                 daily_PAH(n).r_hr_en = r_hr_en;
        %                 daily_PAH(n).week_no = week_no;
    end
end