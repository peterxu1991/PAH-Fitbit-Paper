%%% Fitbit step & Heart rate

clear;clc;close all;

% read the clinic visit date
%% read data

PathName = 'C:\Users\zxu11\OneDrive - Johns Hopkins University\JHU\Projects\Fitbit_IRB_28FEB2018\';
FileName = 'PC_Step_HR_110220.mat';
load([PathName,FileName]);

% pre-allocation of wear pct daily
day_pct = zeros(length(RawFBData),1);
night_pct = zeros(length(RawFBData),1);
sleep_pct = zeros(length(RawFBData),1);

step_avg_hour = zeros(24,2,length(RawFBData));
hr_avg_hour = zeros(24,2,length(RawFBData));
act_avg_hour = zeros(24,2,length(RawFBData));  

% for each individual
for n = 1 : length(RawFBData)
    
    empty_sub = [1,5,8,12,17,18];
    
    if ~isempty(find(empty_sub == n))
        continue;
    else
        
    subi_step = RawFBData(n).step;
    subi_hr = RawFBData(n).hr;
    subi_time = RawFBData(n).time;
    subi_name = RawFBData(n).name;
    subi_act = RawFBData(n).act-1;
    subi_age = RawFBData(n).age;
    
    % use the start and end date to clean up data
    if n == 4 || n == 5 || n ==6 || n==9 ||n == 10 ||n==11
        [subi_hr,subi_step,subi_act,subi_time] = trunc_date(n,subi_hr,subi_step,subi_act,subi_time);
    end
    % remove the repeated rows with same time
    [subi_time,ia, ~] = unique(subi_time);
    subi_hr = subi_hr(ia);
    subi_step = subi_step(ia);
    subi_act = subi_act(ia);
    
    % age-predicted maximum heart rate (refer to Tanaka et al. 2001)
    ap_mhr = 208 - 0.7 * subi_age;
    
    L = 1: 1 : length(subi_time);
    
    % convert time with datenum format to datetime format
    t = datetime(subi_time,'ConvertFrom','datenum');
    
    %% pre-process data
    % remove zero heart rate points
    [rx,~] = find(subi_hr ==0);
    subi_hr(rx) = [];
    subi_step(rx) = [];
    subi_act(rx)=[];
    t(rx) = [];
    
    % keep heart rate data if it is between 20 and age-predicted max heart rate + SD (10bpm) (refer to Avram R et al._2019)
    [rx1,~] = find(subi_hr > ap_mhr + 10 | subi_hr < 20);
    subi_hr(rx1) = [];
    subi_step(rx1) = [];
    subi_act(rx1) = [];
    t(rx1) = [];
    
     %% plot the wear time fraction for each day
        
%      plot_time_fract(t,subi_act,subi_name);
        
    %% plot the missing data distribution
     [week_no, compliance, missing_hours,day_pct_week,night_pct_week,area_step_week,max_hour_step, ...
         core_hour1, core_hour2, step_en, hr_en, r_step_en, r_hr_en] = Fitbit_MissingData2(subi_hr,subi_step,subi_act,t,subi_name,ap_mhr);
   
%% plot the filling data distribution
%             [week_no,mean_step_hour_weekly,std_step_hour_weekly,mean_hr_hour_weekly,std_hr_hour_weekly,...
%         mean_act_hour_weekly,std_act_hour_weekly] = Fitbit_WeekPlot1(subi_hr,subi_step,subi_act,t,subi_name);
%     
%               step_avg_hour(:,:,n) = [mean_step_hour_weekly,std_step_hour_weekly];
%               hr_avg_hour(:,:,n) = [mean_hr_hour_weekly,std_hr_hour_weekly];
%               act_avg_hour(:,:,n) = [mean_act_hour_weekly,std_act_hour_weekly];
    %% Analyze the Density map of activity for each week and
    
    %         [HR_feat, SC_feat, Act_feat]= PAHFitbit_TextureAnalysis(subi_name, week_no);
    %
    %         text_feat(n).HR_feat = HR_feat;
    %         text_feat(n).SC_feat = SC_feat;
    %         text_feat(n).Act_feat = Act_feat;
    %% process the data on the weekly basis at different window size in the temporal space
         [act_abs,act_pct, fa,hr_mu0,hr_sig0,hr_sk0,hr_ku0,hr_ks0,hr_mu1,hr_sig1,...
    hr_sk1,hr_ku1,hr_ks1,sc_mu,sc_sig,sc_sk,sc_ku,sc_ks,slope,intercept,area1,...
    max_hr_week,hr_max_sc,max_step_week,sc_max_hr,avg_top_sc,avg_top_hr,...
    avg_6MWT_sc,avg_6MWT_hr,Pred_SC170,Pred_HR170,weekly_raw] = Fitbit_DataProcess4(subi_hr,subi_step,subi_act,t,subi_name);
    
            %% save the variables to mat file
            daily_Healthy(n).age = subi_age;
            daily_Healthy(n).name = subi_name;
            
            daily_Healthy(n).hr_mu0 = hr_mu0;
            daily_Healthy(n).hr_sig0 = hr_sig0;
            daily_Healthy(n).hr_sk0 = hr_sk0;
            daily_Healthy(n).hr_ku0 = hr_ku0;
            daily_Healthy(n).hr_ks0 = hr_ks0;
            
            daily_Healthy(n).hr_mu1 = hr_mu1;
            daily_Healthy(n).hr_sig1 = hr_sig1;
            daily_Healthy(n).hr_sk1 = hr_sk1;
            daily_Healthy(n).hr_ku1 = hr_ku1;
            daily_Healthy(n).hr_ks1 = hr_ks1;
            
            daily_Healthy(n).sc_mu = sc_mu;
            daily_Healthy(n).sc_sig = sc_sig;
            daily_Healthy(n).sc_sk = sc_sk;
            daily_Healthy(n).sc_ku = sc_ku;
            daily_Healthy(n).sc_ks = sc_ks;
            
            daily_Healthy(n).slope = slope;
            daily_Healthy(n).intercept = intercept;
            daily_Healthy(n).area = area1;   
           
            daily_Healthy(n).act_abs = act_abs;
            daily_Healthy(n).act_pct = act_pct;
          
            daily_Healthy(n).facti = fa;            
            daily_Healthy(n).day_pct = day_pct_week;
            daily_Healthy(n).night_pct = night_pct_week;
            daily_Healthy(n).compliance = compliance;
            daily_Healthy(n).missing_hours = missing_hours;

            daily_Healthy(n).area_step_week = area_step_week;
            daily_Healthy(n).max_hour_step = max_hour_step;
            
            daily_Healthy(n).core_hour1 = core_hour1;
            
            daily_Healthy(n).max_step = max_step_week;
            daily_Healthy(n).hr_max_sc = hr_max_sc;
            daily_Healthy(n).max_hr = max_hr_week;
            daily_Healthy(n).sc_max_hr = sc_max_hr;
            
            daily_Healthy(n).mean_top_sc = avg_top_sc;
            daily_Healthy(n).mean_top_hr = avg_top_hr;
            daily_Healthy(n).MWT_sc = avg_6MWT_sc;
            daily_Healthy(n).MWT_hr = avg_6MWT_hr;
            
            daily_Healthy(n).SC170 = Pred_SC170;
            daily_Healthy(n).HR170 = Pred_HR170;
            daily_Healthy(n).weekly_raw = weekly_raw;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
%             daily_Healthy(n).core_hour2 = core_hour2;
%             
%             daily_Healthy(n).step_en = step_en;
%             daily_Healthy(n).hr_en = hr_en;
%             
%             daily_Healthy(n).r_step_en = r_step_en;
%             daily_Healthy(n).r_hr_en = r_hr_en;
%             daily_Healthy(n).week_no = week_no;
            
    end
end

save 'Daily_Healthy_PC_110220.mat' daily_Healthy;