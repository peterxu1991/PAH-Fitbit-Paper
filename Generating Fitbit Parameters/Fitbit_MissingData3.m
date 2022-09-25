function [week_no,compliance, missing_hours,day_pct_week,night_pct_week,area_step_week,max_hour_step, ...
    core_hour1, core_hour2] = Fitbit_MissingData3(subi_hr,subi_step,subi_act,tt,subi_name,ap_mhr)

%%% Data process
PathName1 = 'C:\Users\zxu11\OneDrive - Johns Hopkins University\JHU\Projects\Fitbit_IRB_28FEB2018\Plots\';

new_folder = [PathName1 'StepDist\' subi_name];
mkdir(new_folder);

new_folder1 = [PathName1 'StepBin\' subi_name];
mkdir(new_folder1);

new_folder2 = [PathName1 'NormStep\' subi_name];
mkdir(new_folder2);
%% calculate the missing time in each week
% no entry
sort_time = sort(unique(tt),'ascend');

[y1,m1,d1] = ymd(sort_time(1));
[yend,mend,dend] = ymd(sort_time(end));

comp_time = datetime(y1,m1,d1,0,0,0): minutes(1):datetime(yend,mend,dend,23,59,0);

[yr,mo,da] = ymd(tt);
[ycomp,~,~] =ymd(comp_time);

no_year = unique(yr);

%% pre-allocation of matrix

% step_fract_week = zeros(24,7,length(no_week));
% wear_fract_week = zeros(24,7,length(no_week));
% hr_fract_week = zeros(24,7,length(no_week));
% act_fract_week = zeros(24,7,length(no_week));

for i = 1 : length(no_year)
    
    ind_year = find(yr == no_year(i));
    ind_year_comp = find(ycomp == no_year(i));
    
    temp_year = tt(ind_year);
    temp_year_comp = comp_time(ind_year_comp);
    
    temp_hr_year = subi_hr(ind_year);
    temp_step_year = subi_step(ind_year);
    %     temp_week = unique(week(temp_year,'weekofyear'));
    temp_week_comp = unique(week(temp_year_comp,'weekofyear'));
    
    for j = 1 : length(temp_week_comp)
        
        ind_week_comp = find(week(temp_year_comp) == temp_week_comp(j));
        ind_week = find(week(temp_year) == temp_week_comp(j));
        
        % make sure the complete week time are correct
        
        if isempty(ind_week)
            
            year(i).compliance_week(j) = NaN;
            year(i).missing_hours_week(j) = 168;
            year(i).max_step_fract(j) = 0;
            year(i).day_pct(j) = 0;
            year(i).night_pct(j) = 0;
            year(i).mean_wear_hour(:,j) = zeros(24,1); year(i).std_wear_hour(:,j) = zeros(24,1);
            year(i).mean_step_hour(:,j) =  zeros(24,1); year(i).std_step_hour(:,j) = zeros(24,1);
            year(i).area_step_week(j) = 0;                
            year(i).max_hour_step(j) = 0;
            year(i).core_hour1(j) = 0;                
            year(i).core_hour2(j) = 0;
            
        else
            
            day_no = weekday(temp_year(ind_week));
            num_day = unique(day_no);
            
            if (i == 1 && j == 1 && range(num_day) < 7) || (i == length(no_year) && j == length(temp_week_comp) && range(num_day) < 7)
                
                year(i).compliance_week(j) = NaN;
                year(i).missing_hours_week(j) = 168;
                year(i).max_step_fract(j) = NaN;
                year(i).day_pct(j) = NaN;
                year(i).night_pct(j) = NaN;
                year(i).mean_wear_hour(:,j) = NaN* ones(24,1); year(i).std_wear_hour(:,j) = NaN* ones(24,1);
                year(i).mean_step_hour(:,j) =  NaN* ones(24,1); year(i).std_step_hour(:,j) = NaN* ones(24,1);
                year(i).area_step_week(j) = NaN;
                year(i).max_hour_step(j) = NaN;
                year(i).core_hour1(j) = NaN;                
                year(i).core_hour2(j) = NaN;
                
                continue;
                
            elseif length(day_no) < 2*24 * 60
                
                year(i).compliance_week(j) = NaN;
                year(i).missing_hours_week(j) = 168;
                year(i).max_step_fract(j) = 0;
                year(i).day_pct(j) = 0;
                year(i).night_pct(j) = 0;
                year(i).mean_wear_hour(:,j) = zeros(24,1); year(i).std_wear_hour(:,j) = zeros(24,1);
                year(i).mean_step_hour(:,j) =  zeros(24,1); year(i).std_step_hour(:,j) = zeros(24,1);
                year(i).area_step_week(j) = 0;
                year(i).max_hour_step(j) = 0;
                year(i).core_hour1(j) = 0;
                year(i).core_hour2(j) = 0;
                               
            else
                
                week_time = temp_year(ind_week);
                week_hr = temp_hr_year(ind_week);
                temp_step = temp_step_year(ind_week);
                week_time_comp = temp_year_comp(ind_week_comp);
                
                % no entry part
                rm_nowear_week = length(week_time_comp) - length(unique(week_time)); %report fraction of no entry
                year(i).compliance_week(j) = (length(week_time_comp) - rm_nowear_week)/length(week_time_comp);
                              
                %         % make new directory
%         PathName1 = 'C:\Users\zxu11\OneDrive - Johns Hopkins University\JHU\Projects\Fitbit_IRB_28FEB2018\Plots\';
%         new_folder = [PathName1 'RawData\' subi_name '\'];
%         mkdir(new_folder)
%         
%         fh1 = figure('units','normalized','position',[0 0 1 1]);
%         subplot(2,1,1)
%         grid minor
%         plot(filled_step,'MarkerEdgeColor','b','LineWidth',2)
%         title('Step Count Raw Signal');
%         xlabel('Minutes','FontSize',25,'FontWeight','bold')
%         ylabel('Step Count','FontSize',25,'FontWeight','bold')
%         
%         subplot(2,1,2)
%         plot(filled_hr,'k','LineWidth',2)
%         title('Heart Rate Raw Signal');
%         xlabel('Minutes','FontSize',25,'FontWeight','bold')
%         ylabel('Heart rate','FontSize',25,'FontWeight','bold')
%         grid on      
%         set(gca,'FontSize',12,'FontWeight','bold');
% 
%         filename1=[new_folder '\'  'week' num2str(k-1 + min(no_week))];
%         print(fh1,'-djpeg',filename1);
%         close(fh1);

                dm = weekday(week_time);
                %                 day_time = week_time(dm);
                
                 %% calculate the daytime and nighttime percentage of wearing during the week
                ht = hms(week_time);
                
                day_ind = find(ht >= 8 & ht < 20);
                night_ind = find(ht >= 20 | ht < 8);
                
                year(i).day_pct(j) = length(day_ind)/length(week_time);
                year(i).night_pct(j) = length(night_ind)/length(week_time);
                                
                %% find the missing hour and compliance of the week
               
                missing_time = week_time(2:end) - week_time(1:end-1);
                missing_gap = hours(max(missing_time));
                
                if week_time(end) == week_time_comp(end) && week_time(1) == week_time_comp(1) 
                                       
                   year(i).missing_hours_week(j) = missing_gap;
                   
                elseif week_time(end) < week_time_comp(end) || week_time(1) > week_time_comp(1) 
                    
                    gap_time_end = hours(week_time_comp(end) - week_time(end));
                    gap_time_start = hours(week_time(1) - week_time_comp(1));
                    
                    max_hours = max([gap_time_end, gap_time_start,missing_gap]);
                    
                    year(i).missing_hours_week(j) =  max_hours;
                                           
                end
            
                %% plot wear time and step count at each hour over each week
                for m = 1 : length(num_day)
                    
                    ind_day1 = find(day_no == num_day(m));
                    
                    step_day = temp_step(ind_day1);
                    total_step_day = sum(step_day);
                    %
                    %             hr_day = hr_M_week(ind_day1);
                    %
                    %             act_day = act_M_week(ind_day1);
                    
                    [hd1,~,~] = hms(week_time(ind_day1));
                    
                    % determine the wear time percentage and average activity level in each hour
                    for n = 1 : 24
                        
                        ind_hr = find(hd1 == n-1);
                        
                        if isempty(ind_hr)
                            
                            wear_fract(n,m) = 0;
                            step_fract(n,m) = NaN;
                            %                     hr_fract(n,i) = NaN;
                            %                     act_fract(n,i) = NaN;
                            
                        else
                            
                            step_fract(n,m) = sum(step_day(ind_hr));
                            %                     hr_fract(j,i) = mean(hr_day(ind_hr));
                            %
                            %                     act1_pct = length(find(act_day(ind_hr) == 1));
                            %                     act2_pct = length(find(act_day(ind_hr) == 2));
                            %                     act3_pct = length(find(act_day(ind_hr) == 3));
                            %                     act4_pct = length(find(act_day(ind_hr) == 4));
                            %
                            %                     act_fract(j,i) = (act1_pct * 1 + act2_pct * 2 + act3_pct * 3 + act4_pct * 4)/length(ind_hr);
                            wear_fract(n,m) = length(ind_hr)/60;
                            
                        end
                    end
                    max_step_fract(m) = max(step_fract(:,m))/total_step_day;
                end
                
                year(i).max_step_fract(j) = nanmean(max_step_fract);
                % calculate mean value for each hour during the week
                mean_wear = nanmean(wear_fract,2); std_wear = nanstd(wear_fract,[],2);
                mean_step = nanmean(step_fract,2); std_step = nanstd(step_fract,[],2);
                %         mean_hr = nanmean(hr_fract,2); std_hr = nanstd(hr_fract,[],2);
                %         mean_act = nanmean(act_fract,2); std_act = nanstd(act_fract,[],2);
                
                year(i).mean_wear_hour(:,j) = mean_wear; year(i).std_wear_hour(:,j) = std_wear;
                year(i).mean_step_hour(:,j) = mean_step; year(i).std_step_hour(:,j) = std_step;
                %             year(i).mean_hr_hour(:,j) = mean_hr; year(i).std_hr_hour(:,j) = std_hr;
                %             year(i).mean_act_hour(:,j) = mean_act; year(i).std_act_hour(:,j) = std_act;
                
                %             step_fract_week(:,num_day,j) = step_fract;
                %             hr_fract_week(:,num_day,j) = hr_fract;
                %             act_fract_week(:,num_day,j) = act_fract;
                
                %% calculate the area covered by mean and std:remove the baseline and normalize the std based on the max std
                rm_mean_step = mean_step;
                norm_std = std_step;
                
                % remove the NaN in the vector
                [nan_step,~] = find(isnan(rm_mean_step));
                rm_mean_step(nan_step) = [];
                norm_std(nan_step) = [];
               area1 =  sum(norm_std);
%                 area2 = polyarea(rm_mean_step,norm_std);
                
                year(i).area_step_week(j) = area1;
                
                year(i).max_hour_step(j) = max(mean_step);
                
                %% find the specific points using two different methods:larger than average of hourly steps on the daytime
                % method 1
                daytime_step = mean_step(8:19);
                
                mean_daytime_step = nanmean(daytime_step);
                
                ind_pks1 = find(mean_step > mean_daytime_step);
                core_step1 = mean_step(ind_pks1);
               
                year(i).core_hour1(j) = length(core_step1);
                                               
                % method 2
                threshold = 0.05;
                baseline = nanmean(mean_step) + threshold * (max(mean_step) - nanmean(mean_step));
                
                [pks_step,locs] = findpeaks(mean_step);
                
                ind_pks2 = find(pks_step > baseline);
                core_step2 = pks_step(ind_pks2);
                ind_locs2 = locs(ind_pks2);
                year(i).core_hour2(j) = length(core_step2);
                %% divide the total step at each hour into different bin size and calculate the fraction during the week
                 temp_step_hourly = step_fract(:);
                ind_nan = find(~isnan(temp_step_hourly));
                all_step_hourly = temp_step_hourly(ind_nan);
                
                bin = 0 : 500 : 5000;
                ct_bin = zeros(11,1);
                
                for s = 1 : 11
                    
                    if s == 11
                        ibin = bin(s);
                        inds = find(all_step_hourly > ibin);
                        ct_bin(s) = length(inds)/length(all_step_hourly);
                    else
                        
                        bin_low = bin(s);
                        bin_high = bin(s+1);
                        inds2 = find(all_step_hourly > bin_low & all_step_hourly < bin_high);
                        ct_bin(s) = length(inds2)/length(all_step_hourly);
                    end
                end
                
                %% plot the fraction of hourly total step distribution
                %          Lname = {'0~500','501~1000','1001~1500','1501~2000','2001~2500','2501~3000','3001~3500','3501~4000','4001~4500','4501~5000','>5000'};
                %         fh1 = figure('units','normalized','position',[0 0 1 1]);
                %         bar(ct_bin);
                %         h = gca;
                %         set(gca,'xticklabel',Lname);
                %         h.XTickLabelRotation = 45;
                %         xlabel('Step Count');
                %         ylabel('Fraction of total hour over the week');
                % %         title(['Total Steps at each hour during the week' num2str(j) 'at year' num2str(i)]);
                %         set(gca,'FontSize',15,'FontWeight','bold');
                %
                %         filename1=[new_folder1 '\'  subi_name '_week' num2str(j) '_year' num2str(i)];
                %         print(fh1,'-djpeg',filename1);
                %         close(fh1);
                
                %% plot the hourly total step distribution
%                         fh1 = figure('units','normalized','position',[0 0 1 1]);
%                         errorbar(mean_step,std_step,'-s','MarkerSize',10,'LineWidth',1.5)
%                         hold on
%                         plot(ind_locs1,core_step1,'o','MarkerSize',12,'MarkerEdgeColor','red','MarkerFaceColor','red','LineWidth',2);
%                         plot(ind_locs2,core_step2,'*','MarkerSize',14,'MarkerEdgeColor','black','MarkerFaceColor','black','LineWidth',2);
%                         xlabel('Hour');
%                         ylabel('Step Count');
%                         title(['Total Steps at each hour during the week' num2str(j) 'at year' num2str(i)]);
%                         legend('Hourly Step','Method 1','Method 2')
%                         set(gca,'FontSize',15,'FontWeight','bold');
%                 
%                         filename1=[new_folder '\'  subi_name '_week' num2str(j) '_year' num2str(i)];
%                         print(fh1,'-djpeg',filename1);
%                         close(fh1);
               %% plot the normalized hourly total step distribution
%                         fh2 = figure('units','normalized','position',[0 0 1 1]);
%                         errorbar(rm_mean_step,norm_std,'-s','MarkerSize',10,'LineWidth',1.5)
%                       
%                         xlabel('Hour');
%                         ylabel('Step Count');
%                         title(['Total Steps at each hour during the week' num2str(j) 'at year' num2str(i)]);
%                         legend('Hourly Step','Method 1','Method 2')
%                         set(gca,'FontSize',15,'FontWeight','bold');
%                 
%                         filename2=[new_folder2 '\'  subi_name '_week' num2str(j) '_year' num2str(i)];
%                         print(fh2,'-djpeg',filename2);
%                         close(fh2);
            end
        end
    end
end


%% plot the missing time component for each week
cell_year = struct2cell(year);

if length(year) > 1
    
    temp_compliance = cell2mat((squeeze(cell_year(1,1,:)))');
    temp_missing_hours = cell2mat((squeeze(cell_year(2,1,:)))');
    
    max_step_fract = cell2mat((squeeze(cell_year(3,1,:)))');
    temp_day_pct_week = cell2mat((squeeze(cell_year(4,1,:)))');
    temp_night_pct_week = cell2mat((squeeze(cell_year(5,1,:)))');
    
    mean_wear_week = cell2mat((squeeze(cell_year(6,1,:)))');
    mean_step_week = cell2mat((squeeze(cell_year(8,1,:)))');
    std_step_week = cell2mat((squeeze(cell_year(9,1,:)))');
    
    temp_area_step_week = cell2mat((squeeze(cell_year(10,1,:)))');
    temp_max_hour_step = cell2mat((squeeze(cell_year(11,1,:)))');
    temp_core_hour1 = cell2mat((squeeze(cell_year(12,1,:)))');
    temp_core_hour2 = cell2mat((squeeze(cell_year(13,1,:)))');
    temp_step_en = cell2mat((squeeze(cell_year(14,1,:)))');
    temp_hr_en = cell2mat((squeeze(cell_year(15,1,:)))');
    temp_r_step = cell2mat((squeeze(cell_year(16,1,:)))');
    temp_r_hr = cell2mat((squeeze(cell_year(17,1,:)))');
                
else
    
    temp_compliance = cell_year{1};
    temp_missing_hours = cell_year{2};
    
    if length(cell_year)>= 3
        
        max_step_fract = cell_year{3};
        temp_day_pct_week = cell_year{4};
        temp_night_pct_week = cell_year{5};
        mean_wear_week = cell_year{6};
        mean_step_week = cell_year{8};
        std_step_week = cell_year{9};
        
        temp_area_step_week = cell_year{10};
        temp_max_hour_step = cell_year{11};
        temp_core_hour1 = cell_year{12};
        temp_core_hour2 = cell_year{13};
        temp_step_en = cell_year{14};
        temp_hr_en = cell_year{15};
        temp_r_step = cell_year{16};
        temp_r_hr = cell_year{17};
    
    end
    
end

% remove Nan in the matrix
ind_num = find(~isnan(temp_compliance));
compliance = (temp_compliance(ind_num))';
missing_hours = (temp_missing_hours(ind_num))';
day_pct_week = (temp_day_pct_week(ind_num))';
night_pct_week = (temp_night_pct_week(ind_num))';
area_step_week = (temp_area_step_week(ind_num))';
max_hour_step = (temp_max_hour_step(ind_num))';
core_hour1 = (temp_core_hour1(ind_num))';
core_hour2 = (temp_core_hour2(ind_num))';

% determine the true week no without excluding the missing gaps
[~,~,week_no] = intersect(compliance,temp_compliance,'stable');
%% plot the area covered by mean and upper std at the hour each week
%         fh1 = figure('units','normalized','position',[0 0 1 1]);
% 
%         plot(area_step_week,'d','MarkerSize',12,'LineWidth',2);         
%         grid minor
%         xlabel('Week No');
%         ylabel('Area');
%         title([subi_name],'Interpreter','none');
%         set(gca,'FontSize',15,'FontWeight','bold');
% 
%         filename1=[PathName1 'AreaStep\' subi_name '_AreaStep'];
%         print(fh1,'-djpeg',filename1);
%         close(fh1);
%% plot the fraction of max step at the hour over the total step each day
%         fh1 = figure('units','normalized','position',[0 0 1 1]);
%
%         plot(max_step_fract,'d','MarkerSize',12,'LineWidth',2);
%
%         ylim([0,1]);
%         set(gca,'YTick',0:0.1:1)
%         grid minor
%         xlabel('Week No');
%         ylabel('Fraction of Max Hourly Step during the day');
%         title([subi_name],'Interpreter','none');
%         set(gca,'FontSize',15,'FontWeight','bold');
%
%         filename1=[PathName1 'MaxStep\' subi_name '_MaxStep'];
%         print(fh1,'-djpeg',filename1);
%         close(fh1);
%% plot the missing hours vs compliance
%         fh1 = figure('units','normalized','position',[0 0 1 1]);
%         [indx,indy] = find(missing_hours == 0);
%         missing_hours(indx,indy) = 1;
%         log_miss = log10(missing_hours);
% %         data1 = [compliance(1:end-1),log_miss(1:end-1)];
% %         data2 = [compliance(2:end),log_miss(2:end)];
%         x = compliance;
%         y = log_miss;
%         dx = range(x)/100;
%         dy = range(y)/100;
%        
%         xno = num2str(week_no');
%         pt_name = cellstr(xno);
%                  
% %         arrow3(data1,data2);
%         scatter(x,y,100,'filled');
%         text(x+dx, y+dy, pt_name,'FontSize',15);
%         xlim([0,1]);
%         set(gca,'XTick',0:0.1:1);
%         ylim([0,log10(168)]);
%        set(gca,'YTick',0:0.25:log10(168)); 
%        
%         grid minor
%         xlabel('Compliance');
%         ylabel('Longest block without wearing (log)');
%         title([subi_name ' weekly data missing distribution'],'Interpreter','none');
%         set(gca,'FontSize',15,'FontWeight','bold');
% 
%         filename1=[PathName1 'MissingData\MissingHours\' subi_name '_missing_hours'];
%         print(fh1,'-djpeg',filename1);
%         close(fh1);
        
        %% plot the missing hours vs compliance
%         fh2 = figure('units','normalized','position',[0 0 1 1]);
% 
%         x2 = compliance;
%         y2 = missing_hours/168;
%         dx2 = range(x2)/100;
%         dy2 = range(y2)/100;
%                  
%         scatter(x2,y2,100,'filled');
%         text(x2+dx2, y2+dy2, pt_name,'FontSize',15);
%         xlim([0,1]);
%         set(gca,'XTick',0:0.1:1);
%         ylim([0,1]);
%        set(gca,'YTick',0:0.1:1); 
% %         ylim([0,25]);
% %         
%         grid minor
%         xlabel('Compliance');
%         ylabel('Longest block without wearing time/total hours of complete week');
%         title([subi_name ' weekly data missing distribution'],'Interpreter','none');
%         set(gca,'FontSize',15,'FontWeight','bold');
% 
%         filename2=[PathName1 'MissingData\MissingHours1\' subi_name '_missing_hours'];
%         print(fh2,'-djpeg',filename2);
%         close(fh2);
%% plot the hourly total step distribution
%         fh1 = figure('units','normalized','position',[0 0 1 1]);
%         errorbar(mean_step_week,std_step_week,'-s','MarkerSize',10,...
%             'MarkerEdgeColor','red','MarkerFaceColor','red','LineWidth',2);
%         xlabel('Hour');
%         ylabel('Step Count');
%         title(['Total Steps at each hour during the week']);
%         set(gca,'FontSize',15,'FontWeight','bold');
%
%         filename1=[PathName1 'StepDist\' subi_name '_week'];
%         print(fh1,'-djpeg',filename1);
%         close(fh1);
%% plot the wear time heat map based on the hourly percentage of wear time over each week
%         fh1 = figure('units','normalized','position',[0 0 1 1]);
%         % tb1 = table(wear_pct_hr);
%         hmap = heatmap(mean_wear_week);
%         hmap.XLabel = 'Week';
%         hmap.YLabel = 'Hour';
%         hmap.Colormap = parula;
%
%         filename1=[PathName1 'WearTime_Week\' subi_name 'WearTime_map'];
%         print(fh1,'-djpeg',filename1);
%         close(fh1);
end
