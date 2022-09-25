function plot_time_fract(t, subi_act, subi_name)

% To distinguish the activity level = 0 with the N/A cell, add 1 to the activity

[maxY,maxM,maxD] = ymd(max(t));
[minY,minM,minD] = ymd(min(t));

max_date = datetime(maxY,maxM,maxD);
min_date = datetime(minY,minM,minD);

length_day = daysact(min_date,max_date);

wear_pct_fract = zeros(24,length_day+1);
act_fract = zeros(24,length_day+1); 

[yt,mt,dt] = ymd(t);

for i = 1 : length_day+1
    
    dayi = min_date + days(i)-days(1);
    [yi,mi,di] = ymd(dayi);
    
    ind_day = find(yt == yi & mt == mi & dt == di);
    [hd,md,~] = hms(t(ind_day));
    act_day = subi_act(ind_day);
    
    no_h = unique(hd);
    
%     % Determine the activity level in each minute
%     idxt = hd*60 + md + 1;
%     
%     act_fract(idxt,i) = act_day;

    % determine the wear time percentage and average activity level in each hour
    for j = 1 : 24
        
        ind_hr = find(hd == j-1);
              
        if isempty(ind_hr)
            wear_pct_fract(j,i) = NaN;
            act_fract(j,i) = NaN;
                        
        else
            
            act1_pct = length(find(act_day(ind_hr) == 1));
            act2_pct = length(find(act_day(ind_hr) == 2));
            act3_pct = length(find(act_day(ind_hr) == 3));
            act4_pct = length(find(act_day(ind_hr) == 4));
            
            act_fract(j,i) = (act1_pct * 1 + act2_pct * 2 + act3_pct * 3 + act4_pct * 4)/length(ind_hr);
             wear_pct_fract(j,i) = length(ind_hr)/60;
        end
        
        
        
    end
end

%% plot the wear time heat map based on the hourly percentage of wear time
fh1 = figure('units','normalized','position',[0 0 1 1]);
% tb1 = table(wear_pct_hr);
hmap = heatmap(wear_pct_fract);
hmap.XLabel = 'Day';
hmap.YLabel = 'Hour';
hmap.Colormap = parula;
PathName1 = 'C:\Users\zxu11\OneDrive - Johns Hopkins University\JHU\Projects\Fitbit_IRB_28FEB2018\Plots\WearTime\';

filename1=[PathName1 subi_name 'WearTime_map'];
print(fh1,'-djpeg',filename1);
close(fh1);

%% plot the activity level heat map at each minute
% fh2 = figure('units','normalized','position',[0 0 1 1]);
% hmap = heatmap(act_fract);
% caxis(hmap,[0,3]);
% hmap.XLabel = 'Day';
% hmap.YLabel = 'Activity Level';
% hmap.Colormap = parula;
% 
% PathName2 = 'C:\Users\zxu11\OneDrive - Johns Hopkins University\JHU\Projects\Fitbit_IRB_28FEB2018\COPD_Plots\ActLevel\';
% 
% filename2=[PathName2 subi_name 'ActLevel_map'];
% print(fh2,'-djpeg',filename2);
% close(fh2);
