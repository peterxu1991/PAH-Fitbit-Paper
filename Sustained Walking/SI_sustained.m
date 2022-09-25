%% Clear workspace

clear;clc;close all;

%% Organize data for healthy, PC, and COPD patients

load('Daily_Healthy_PC_021521.mat'); % This structure contains data for healthy, PC, and COPD patients.
load('daily_PAH_082520.mat'); % This structure contains data for PAH patients.

% List indices for healthy, PC, and COPD patients in structure
healthy_indices = [2; 3; 7; 9; 13; 19];
PC_indices = [4; 6];
COPD_indices = [10, 11, 14:18, 20:32];

% Create new structures for each patient population.
weekly_raw_healthy = struct([]);
weekly_raw_PC = struct([]);
weekly_raw_COPD = struct([]);

% Compile weekly_raw data for each healthy subject.
for i = 1:length(healthy_indices)
    weekly_raw_healthy(i).name = daily_Healthy(healthy_indices(i)).name;
    weekly_raw_healthy(i).weekly_raw = struct([]);
    for j = 1:length(daily_Healthy(healthy_indices(i)).weekly_raw)
        if isempty(daily_Healthy(healthy_indices(i)).weekly_raw(j).YearNo)
        else
            weekly_raw_healthy(i).weekly_raw = [weekly_raw_healthy(i).weekly_raw daily_Healthy(healthy_indices(i)).weekly_raw(j)];
        end
    end
end

% Compile weekly_raw data for each PC patient.
for i = 1:length(PC_indices)
    weekly_raw_PC(i).name = daily_Healthy(PC_indices(i)).name;
    weekly_raw_PC(i).weekly_raw = struct([]);
    for j = 1:length(daily_Healthy(PC_indices(i)).weekly_raw)
        if isempty(daily_Healthy(PC_indices(i)).weekly_raw(j).YearNo)
        else
            weekly_raw_PC(i).weekly_raw = [weekly_raw_PC(i).weekly_raw daily_Healthy(PC_indices(i)).weekly_raw(j)];
        end
    end
end

% Compile weekly_raw data for each COPD patient.
for i = 1:length(COPD_indices)
    weekly_raw_COPD(i).name = daily_Healthy(COPD_indices(i)).name;
    weekly_raw_COPD(i).weekly_raw = struct([]);
    for j = 1:length(daily_Healthy(COPD_indices(i)).weekly_raw)
        if isempty(daily_Healthy(COPD_indices(i)).weekly_raw(j).YearNo)
        else
            weekly_raw_COPD(i).weekly_raw = [weekly_raw_COPD(i).weekly_raw daily_Healthy(COPD_indices(i)).weekly_raw(j)];
        end
    end
end

%%% General Parameters
SCthresh = 60; % Count sustained min of walking for 60+ SPM
histdatathresh = 5; % Need at least 5 data points (of sustained walking) to include in plots
min_dur = 2; % minimum # of minutes to be counted as "sustained" walking


%% Find # of minutes of sustained walking for healthy subjects

this_weekly_raw = weekly_raw_healthy;

% Gather # of minutes of sustained walking & avg SPM on week-by-week basis
for i = 1:length(this_weekly_raw)
    this_weekly_raw(i).sustained = struct([]);
    
    for j = 1:length(this_weekly_raw(i).weekly_raw)
        this_weekly_raw(i).sustained(j).WeekNo = this_weekly_raw(i).weekly_raw(j).WeekNo;
        this_weekly_raw(i).sustained(j).YearNo = this_weekly_raw(i).weekly_raw(j).YearNo;
        this_weekly_raw(i).sustained(j).duration = [];
        this_weekly_raw(i).sustained(j).avgSPM = [];
        this_weekly_raw(i).sustained(j).rawdata = struct([]);
        
        this_week = this_weekly_raw(i).weekly_raw(j);     
        this_time = this_week.all_time;

        time_start = datevec(datenum(this_time(1))); % Time at start of week
        time_end = datevec(datenum(this_time)); % Time at end of week
        elapsed_time = etime(time_end,time_start); % in seconds

        this_SC = this_week.all_SC;
        this_HR = this_week.all_HR; % Unused for now (2/25/2021).
        % [h,m,s] = hms(this_time);
        
        suswalk_idx = [];
        
        for k = 1:length(this_SC)-1
            if this_SC(k) > SCthresh
                streak = 1;
                while k+streak < length(this_SC)+1
                    if this_SC(k+streak) > SCthresh
                        if elapsed_time(k+streak) - elapsed_time(k) == streak*60
                            streak = streak + 1;
                        else
                            break
                        end
                    else
                        break
                    end
                end
                
                if k+streak > length(this_SC) % For edge case where sustained walking occurs at the very end of vector
                    streak = streak - 1;
                end
                
                if streak == 1
                    continue
                else
                    suswalk_idx = [suswalk_idx; k];
                    this_weekly_raw(i).sustained(j).duration = [this_weekly_raw(i).sustained(j).duration; streak];
                    this_weekly_raw(i).sustained(j).avgSPM = [this_weekly_raw(i).sustained(j).avgSPM; mean(this_SC(k:k+streak-1))];
                end
            end
        end
        
        % Remove "streaks" of sustained walking that are within the same ambulation.
        final_duration = [];
        final_avgSPM = [];
        final_suswalk_idx = [];
        idx = 1;
        while idx < length(this_weekly_raw(i).sustained(j).duration)+1
            this_duration = this_weekly_raw(i).sustained(j).duration(idx);
            final_duration = [final_duration; this_duration];
            final_avgSPM = [final_avgSPM; this_weekly_raw(i).sustained(j).avgSPM(idx+this_duration-2)];
            final_suswalk_idx = [final_suswalk_idx; suswalk_idx(idx)];
            idx = idx + this_weekly_raw(i).sustained(j).duration(idx) - 1;
        end
        this_weekly_raw(i).sustained(j).final_duration = final_duration; 
        this_weekly_raw(i).sustained(j).final_avgSPM = final_avgSPM;
        this_weekly_raw(i).sustained(j).compileddata = horzcat(final_duration, final_avgSPM);
        
        % Calculate 1/e value of weekly duration histogram
        dur_data = this_weekly_raw(i).sustained(j).compileddata(:,1);
        dur_data_mod = [];
        for k = 1:length(dur_data)
            if dur_data(k) < min_dur
            else
                dur_data_mod = [dur_data_mod; dur_data(k)];
            end
        end
        
        if length(unique(dur_data_mod)) > 1
            this_weekly_raw(i).sustained(j).dist_var = fitdist(dur_data_mod-min_dur,'exponential');
            this_weekly_raw(i).sustained(j).expdist_mu = this_weekly_raw(i).sustained(j).dist_var(1).mu+min_dur;
        end
        
        for k = 1:length(final_suswalk_idx)
            % Store raw data
            start_idx = final_suswalk_idx(k);
            end_idx = start_idx + final_duration(k) - 1;
            this_weekly_raw(i).sustained(j).rawdata(k).rawSPM = this_week.all_SC(start_idx:end_idx);
            this_weekly_raw(i).sustained(j).rawdata(k).rawHR = this_week.all_HR(start_idx:end_idx);
        end
    end
    
    temp1 = []; temp2 = [];
    for j = 1:length(this_weekly_raw(i).sustained)
        temp1 = [temp1; this_weekly_raw(i).sustained(j).expdist_mu];
        temp2 = [temp2; length(this_weekly_raw(i).sustained(j).final_duration)];
    end
    this_weekly_raw(i).avg_end_expdist_mu = mean(temp1);
    this_weekly_raw(i).avg_freq = mean(temp2);
end

% Compile # of minutes of sustained walking & avg SPM across all weeks
for i = 1:length(this_weekly_raw)
    this_weekly_raw(i).compileddata = [];
    temp_duration = []; temp_avgSPM = [];
    
    for j = 1:length(this_weekly_raw(i).sustained)
        temp_duration = [temp_duration; this_weekly_raw(i).sustained(j).final_duration];
        temp_avgSPM = [temp_avgSPM; this_weekly_raw(i).sustained(j).final_avgSPM];
    end
    
    this_weekly_raw(i).compileddata = horzcat(temp_duration, temp_avgSPM);
    
end

% Find total number of sustained steps for each week
for i = 1:length(this_weekly_raw)
    for j = 1:length(this_weekly_raw(i).sustained)
        temp = this_weekly_raw(i).sustained(j).final_duration .* this_weekly_raw(i).sustained(j).final_avgSPM;
        this_weekly_raw(i).sustained(j).totalsussteps = sum(temp);
        this_weekly_raw(i).sustained(j).totalsteps = sum(this_weekly_raw(i).weekly_raw(j).all_SC);
        this_weekly_raw(i).sustained(j).susovertotal = this_weekly_raw(i).sustained(j).totalsussteps/this_weekly_raw(i).sustained(j).totalsteps;
    end
end

weekly_raw_healthy = this_weekly_raw;


%% Statistical Analysis of 3D Histogram/Heat Map

this_weekly_raw = weekly_raw_healthy;

for i = 1:length(this_weekly_raw)
    % this_weekly_raw(i).totalsteps_week = mean(this_weekly_raw(i).norm_SC(1).weekly_norm_meanSC)*7;
    
    temp1 = []; temp2 = []; temp3 = [];
    for j = 1:length(this_weekly_raw(i).sustained)
        temp1 = [temp1; this_weekly_raw(i).sustained(j).totalsussteps];
        temp2 = [temp2; this_weekly_raw(i).sustained(j).totalsteps];
        temp3 = [temp3; this_weekly_raw(i).sustained(j).susovertotal];
    end
    this_weekly_raw(i).avg_totalsussteps = mean(temp1);
    this_weekly_raw(i).avg_totalsteps = mean(temp2);
    this_weekly_raw(i).avg_susovertotal = mean(temp3);
end

weekly_raw_healthy = this_weekly_raw;


%% Graph heat map/3D histogram of sustained walking data for healthy subjects

this_weekly_raw = weekly_raw_healthy;

for i = 1:length(this_weekly_raw)
    if size(this_weekly_raw(i).compileddata,1) > histdatathresh
        % Specify centers of histogram bins
        bins_duration_max = max(this_weekly_raw(i).compileddata(:,1));
        % bins_duration_max = prctile(this_weekly_raw(i).compileddata(:,1),95);
        bins_avgSPM_max = max(this_weekly_raw(i).compileddata(:,2));
        bins_duration = 2:1:bins_duration_max;
        bins_duration_2D = 2:1:(bins_duration_max+1);
        bins_duration_2D = bins_duration_2D - 0.5;
        bins_avgSPM = 60:10:bins_avgSPM_max;
        
        [this_weekly_raw(i).N, this_weekly_raw(i).C] = hist3(this_weekly_raw(i).compileddata,'Ctrs',{bins_duration bins_avgSPM});

%         %%% Plot 3D histogram (sustained walking duration & avg SPM)
%         figure(i)
%         hold on
%         [this_weekly_raw(i).histcount,this_weekly_raw(i).histcenter] = hist3(this_weekly_raw(i).compileddata,'Ctrs',{bins_duration bins_avgSPM});
%         hist3(this_weekly_raw(i).compileddata,'Ctrs',{bins_duration bins_avgSPM},'CdataMode','auto')
%         title(sprintf('Histogram: Sustained Walking Duration & Avg SPM for %s',string(this_weekly_raw(i).name)))
%         xlabel('Duration (min)')
%         ylabel('SPM')
%         % colorbar
%         view(3)
%         TextLocation(['Avg # of sustained steps/week = ' num2str(sprintf('%d',round(this_weekly_raw(i).avg_totalsteps)))],...
%             'Location','best');
%         hold off
% 
%         %%% Plot heat map (sustained walking duration & avg SPM)
%         figure(i+length(this_weekly_raw))
%         hold on
%         hist3(this_weekly_raw(i).compileddata,'Ctrs',{bins_duration bins_avgSPM},'CdataMode','auto')
%         title(sprintf('Heat Map: Sustained Walking Duration & Avg SPM for %s',string(this_weekly_raw(i).name)))
%         xlabel('Duration (min)')
%         ylabel('SPM')
%         colorbar
%         hold off
%         
%         %%% Plot 2D duration histogram (disregard avg SPM)
%         figure(i+2*length(this_weekly_raw))
%         hold on
%         histogram(this_weekly_raw(i).compileddata(:,1),bins_duration_2D)
%         % Statistics
%         dur_median = median(this_weekly_raw(i).compileddata(:,1));
%         dur_sig = std(this_weekly_raw(i).compileddata(:,1));
%         dur_sk = skewness(this_weekly_raw(i).compileddata(:,1));
%         [~,~,dur_ks,~] = kstest(this_weekly_raw(i).compileddata(:,1),'CDF',fitdist(this_weekly_raw(i).compileddata(:,1),'lognormal'));
%         TextLocation(['Median = ' num2str(sprintf('%.3g',dur_median)), ', SD = ' num2str(sprintf('%.3g',dur_sig)),...
%             newline, 'Skewness = ' num2str(sprintf('%.3g',dur_sk)), ', KS-statistics = ' num2str(sprintf('%.3g',dur_ks))],...
%             'Location','northeast');
%         title(sprintf('Sustained Walking Duration for %s',string(this_weekly_raw(i).name)))
%         xlabel('Duration (min)')
%         ylabel('Frequency')   
%         hold off
        
        %%% Plot modified 2D duration histogram: require minimum # of sustained minutes
        dur_data = this_weekly_raw(i).compileddata(:,1);
        dur_data_mod = [];
        for j = 1:length(dur_data)
            if dur_data(j) < min_dur
            else
                dur_data_mod = [dur_data_mod; dur_data(j)];
            end
        end
        
        if length(unique(dur_data_mod)) > 1 % arbitrary cutoff
            % ~isempty(dur_data_mod)
            figure(i+3*length(this_weekly_raw))
            hold on
            histfit(dur_data_mod-min_dur,length(bins_duration_2D),'exponential')
            % histogram(dur_data_mod,length(bins_duration_2D))
            this_weekly_raw(i).end_dist_var = fitdist(dur_data_mod-min_dur,'exponential');
            
            % Statistics
            dur_median = median(dur_data_mod);
            dur_sig = std(dur_data_mod);
            dur_sk = skewness(dur_data_mod);
            if length(unique(dur_data_mod)) < 2
                dur_ks = NaN;
            else
                [~,~,dur_ks,~] = kstest(dur_data_mod,'CDF',fitdist(dur_data_mod,'lognormal'));
            end
            % 95th percentile of exponential distribution (PDF)
            ninetyfive = -this_weekly_raw(i).end_dist_var(1).mu * log(0.05);
            
            TextLocation(['Median = ' num2str(sprintf('%.3g',dur_median)), ', SD = ' num2str(sprintf('%.3g',dur_sig)),...
                newline, 'Skewness = ' num2str(sprintf('%.3g',dur_sk)), ', KS-statistics = ' num2str(sprintf('%.3g',dur_ks)),...
                newline, newline, '1/e value = ' num2str(sprintf('%.3g',this_weekly_raw(i).end_dist_var(1).mu)),...
                ', 95th percentile = ' num2str(sprintf('%.3g',ninetyfive))],...
                'Location','northeast');
            title(sprintf('Sustained Walking Duration for %s',string(this_weekly_raw(i).name)))
            xlabel('Duration (min)')
            ylabel('Frequency')   
            hold off
        end
        
%         %%% Plot 2D avg SPM histogram (disregard duration)
%         %%% THIS DIDN'T REALLY WORK THOUGH
%         figure(i+4*length(this_weekly_raw))
%         hold on
%         this_data = this_weekly_raw(i).compileddata(:,2);
%         histogram(this_data)
%         % Statistics
%         med = median(this_data);
%         sig = std(this_data);
%         sk = skewness(this_data);
%         [~,~,ks,~] = kstest(this_data,'CDF',fitdist(this_data,'exponential'));
%         TextLocation(['Median = ' num2str(sprintf('%.3g',med)), ', SD = ' num2str(sprintf('%.3g',sig)),...
%             newline, 'Skewness = ' num2str(sprintf('%.3g',sk)), ', KS-statistics = ' num2str(sprintf('%.3g',ks))],...
%             'Location','northeast');
%         title(sprintf('Avg SPM during Sustained Walking for %s',string(this_weekly_raw(i).name)))
%         xlabel('SPM')
%         ylabel('Frequency')   
%         hold off

    end

end

weekly_raw_healthy = this_weekly_raw;


%% Graph scatter plot of endurance vs ambulation frequency for healthy subjects

this_weekly_raw = weekly_raw_healthy;

for i = 1:length(this_weekly_raw)
    xdata = []; ydata = [];
    for j = 1:length(this_weekly_raw(i).sustained)
        if ~isempty(this_weekly_raw(i).sustained(j).expdist_mu)
            xdata = [xdata; length(this_weekly_raw(i).sustained(j).final_duration)];
            ydata = [ydata; this_weekly_raw(i).sustained(j).expdist_mu];
        end
    end
    figure(i)
    hold on
    scatter(xdata,ydata)
    [p,S] = polyfit(xdata, ydata, 1);
    y = @(x) (x*p(1) + p(2));
    fplot(y,'-k','LineWidth',1.5)
    title(sprintf('Endurance (1/e) vs Ambulation Frequency for %s',string(this_weekly_raw(i).name)))
    xlabel('Ambulation Frequency')
    ylabel('Endurance (1/e of exp fit of duration histogram) (min)')
    TextLocation(['y = ' num2str(sprintf('%.3g',p(1))), '*x + ' num2str(sprintf('%.3g',p(2))),...
        newline,newline,'Norm of residuals = ' num2str(sprintf('%.3g',S.normr))],...
        'Location','northeast');
    hold off
end


%% SC Histogram: "Intensity" Parameter for Healthy Subjects (exponential fit)

% this_weekly_raw = weekly_raw_healthy;
% 
% for i = 1:length(this_weekly_raw)
%     this_weekly_raw(i).SChist = struct([]);
%     
%     for j = 1:length(this_weekly_raw(i).weekly_raw)
%         this_weekly_raw(i).SChist(j).WeekNo = this_weekly_raw(i).weekly_raw(j).WeekNo;
%         this_weekly_raw(i).SChist(j).YearNo = this_weekly_raw(i).weekly_raw(j).YearNo;
%         this_weekly_raw(i).SChist(j).dist_var = fitdist(this_weekly_raw(i).weekly_raw(j).SC,'exponential');
%         this_weekly_raw(i).SChist(j).expdist_mu = this_weekly_raw(i).SChist(j).dist_var(1).mu;
%     end
%     
%     temp1 = []; temp2 = [];
%     for j = 1:length(this_weekly_raw(i).SChist)
%         temp1 = [temp1; this_weekly_raw(i).SChist(j).expdist_mu];
%         temp2 = [temp2; this_weekly_raw(i).weekly_raw(j).SC];
%     end
%     this_weekly_raw(i).avg_int_expdist_mu = mean(temp1);
%     this_weekly_raw(i).std_int_expdist_mu = std(temp1);
%     
% %     figure(i)
% %     hold on
% %     histfit(temp2,max(temp2),'exponential')
% %     % Statistics
% %     med = median(temp2);
% %     sig = std(temp2);
% %     sk = skewness(temp2);
% %     if length(unique(temp2)) < 2
% %         ks = NaN;
% %     else
% %         [~,~,ks,~] = kstest(temp2,'CDF',fitdist(temp2,'exponential'));
% %     end
% %                
% %     TextLocation(['Median = ' num2str(sprintf('%.3g',med)), ', SD = ' num2str(sprintf('%.3g',sig)),...
% %         newline, 'Skewness = ' num2str(sprintf('%.3g',sk)), ', KS-statistics = ' num2str(sprintf('%.3g',ks)),...
% %         newline, newline, 'Avg 1/e value = ' num2str(sprintf('%.3g',this_weekly_raw(i).avg_expdist_mu))],...
% %         'Location','northeast');
% %     title(sprintf('Consolidated SC Histogram for %s',string(this_weekly_raw(i).name)))
% %     xlabel('SPM')
% %     ylabel('Frequency')   
% %     hold off
%     
% end
% 
% weekly_raw_healthy = this_weekly_raw;


%% SC Histogram: "Intensity" Parameter for healthy subjects (lognormal fit)

this_weekly_raw = weekly_raw_healthy;

for i = 1:length(this_weekly_raw)
    this_weekly_raw(i).SChist = struct([]);
    
    for j = 1:length(this_weekly_raw(i).weekly_raw)
        this_weekly_raw(i).SChist(j).WeekNo = this_weekly_raw(i).weekly_raw(j).WeekNo;
        this_weekly_raw(i).SChist(j).YearNo = this_weekly_raw(i).weekly_raw(j).YearNo;
        this_weekly_raw(i).SChist(j).dist_var = fitdist(this_weekly_raw(i).weekly_raw(j).SC,'lognormal');
        this_weekly_raw(i).SChist(j).lndist_mu = this_weekly_raw(i).SChist(j).dist_var(1).mu;
        % this_weekly_raw(i).SChist(j).lndist_sigma = this_weekly_raw(i).SChist(j).dist_var(1).sigma;
    end
    
    temp1 = []; temp2 = []; temp3 = [];
    for j = 1:length(this_weekly_raw(i).SChist)
        temp1 = [temp1; this_weekly_raw(i).SChist(j).lndist_mu];
        temp2 = [temp2; this_weekly_raw(i).weekly_raw(j).SC];
        % temp3 = [temp3; this_weekly_raw(i).SChist(j).lndist_sigma];
    end
    this_weekly_raw(i).avg_int_lndist_mu = mean(temp1);
    this_weekly_raw(i).std_int_lndist_mu = std(temp1);
    % this_weekly_raw(i).avg_int_lndist_sigma = mean(temp3);
    
    figure(i)
    hold on
    histfit(temp2,max(temp2),'lognormal')
    pd = fitdist(temp2,'lognormal');
    % Statistics
    med = median(temp2);
    sig = std(temp2);
    this_weekly_raw(i).avg_int_sigma = sig; % Save this into structure
    sk = skewness(temp2);
    if length(unique(temp2)) < 2
        ks = NaN;
    else
        [~,~,ks,~] = kstest(temp2,'CDF',fitdist(temp2,'lognormal'));
    end
               
    TextLocation(['Median = ' num2str(sprintf('%.3g',med)), ', SD = ' num2str(sprintf('%.3g',sig)),...
        newline, 'Skewness = ' num2str(sprintf('%.3g',sk)), ', KS-statistics = ' num2str(sprintf('%.3g',ks))],...
        'Location','northeast');
    title(sprintf('Consolidated SC Histogram for %s',string(this_weekly_raw(i).name)))
    xlabel('SPM')
    ylabel('Frequency')   
    hold off
    
end

weekly_raw_healthy = this_weekly_raw;


%% Graph scatter plot of intensity vs ambulation frequency for healthy subjects

this_weekly_raw = weekly_raw_healthy;

for i = 1:length(this_weekly_raw)
    xdata = []; ydata = [];
    for j = 1:length(this_weekly_raw(i).SChist)
        if ~isempty(this_weekly_raw(i).sustained(j).expdist_mu)
            xdata = [xdata; this_weekly_raw(i).sustained(j).expdist_mu];
            ydata = [ydata; this_weekly_raw(i).SChist(j).expdist_mu];
        end
    end
    figure(i)
    hold on
    scatter(xdata,ydata)
    [p,S] = polyfit(xdata, ydata, 1);
    y = @(x) (x*p(1) + p(2));
    fplot(y,'-k','LineWidth',1.5)
    title(sprintf('Intensity (1/e) vs Endurance (1/e) for %s',string(this_weekly_raw(i).name)))
    xlabel('Endurance (1/e of exp fit of SC histogram) (min)')
    ylabel('Intensity (1/e of exp fit of SC histogram) (SPM)')
    TextLocation(['y = ' num2str(sprintf('%.3g',p(1))), '*x + ' num2str(sprintf('%.3g',p(2))),...
        newline,newline,'Norm of residuals = ' num2str(sprintf('%.3g',S.normr))],...
        'Location','northeast');
    hold off
end


%% Graph scatter plot of intensity vs ambulation frequency for healthy subjects

this_weekly_raw = weekly_raw_healthy;

for i = 1:length(this_weekly_raw)
    xdata = []; ydata = [];
    for j = 1:length(this_weekly_raw(i).SChist)
        if ~isempty(this_weekly_raw(i).SChist(j).expdist_mu)
            xdata = [xdata; length(this_weekly_raw(i).sustained(j).final_duration)];
            ydata = [ydata; this_weekly_raw(i).SChist(j).expdist_mu];
        end
    end
    figure(i)
    hold on
    scatter(xdata,ydata)
    [p,S] = polyfit(xdata, ydata, 1);
    y = @(x) (x*p(1) + p(2));
    fplot(y,'-k','LineWidth',1.5)
    title(sprintf('Intensity (1/e) vs Ambulation Frequency for %s',string(this_weekly_raw(i).name)))
    xlabel('Ambulation Frequency')
    ylabel('Intensity (1/e of exp fit of SC histogram) (SPM)')
    TextLocation(['y = ' num2str(sprintf('%.3g',p(1))), '*x + ' num2str(sprintf('%.3g',p(2))),...
        newline,newline,'Norm of residuals = ' num2str(sprintf('%.3g',S.normr))],...
        'Location','northeast');
    hold off
end


%% Find avg steps per day (using data between 6am - 10pm) for Healthy Subjects

% this_weekly_raw = weekly_raw_healthy;
% 
% % If patient wore the Fitbit the entire time between 6am-10pm, how many data
% % points is that?
% maxwear = 60*16; % minutes
% 
% for i = 1:length(this_weekly_raw)
%     this_weekly_raw(i).norm_SC = struct([]);
%     week = []; year = [];
%     ydata = [];
%     
%     weekly_meanSC = [];
%     weekly_meanwear = [];
%     
%     for j = 1:length(this_weekly_raw(i).weekly_raw)
%         thisweek_days = [];
%         if isempty(this_weekly_raw(i).weekly_raw(j).all_time)
%         else
%             week = [week; this_weekly_raw(i).weekly_raw(j).WeekNo];
%             year = [year; this_weekly_raw(i).weekly_raw(j).YearNo];
%             for k = 1:length(this_weekly_raw(i).weekly_raw(j).all_time)
%                 thisweek_days = [thisweek_days; day(this_weekly_raw(i).weekly_raw(j).all_time(k))];
%             end
%             thisweek_uniquedays = unique(thisweek_days);
%             thisweek_SC = zeros(length(thisweek_uniquedays),3);
%             for k = 1:length(thisweek_uniquedays)
%                 fulldayentries = find(day(this_weekly_raw(i).weekly_raw(j).all_time) == ...
%                     thisweek_uniquedays(k));
%                 
%                 daytimeentries = [];
%                 daytimeSC = [];
%                 
%                 for l = fulldayentries(1):fulldayentries(length(fulldayentries))
%                     if hour(this_weekly_raw(i).weekly_raw(j).all_time(l)) > 5 &&... 
%                             hour(this_weekly_raw(i).weekly_raw(j).all_time(l)) < 22
%                         daytimeentries = [daytimeentries; l];
%                     else
%                     end
%                 end
%                 for l = 1:length(daytimeentries)
%                     daytimeSC = [daytimeSC; this_weekly_raw(i).weekly_raw(j).all_SC(daytimeentries(l))];
%                 end
%                 total_daytimeSC = sum(daytimeSC);
%                 total_daytime_weartime = length(daytimeSC);
%                 thisweek_SC(k,1) = thisweek_uniquedays(k);
%                 thisweek_SC(k,2) = total_daytimeSC;
%                 thisweek_SC(k,3) = total_daytime_weartime;
%             end
%             weekly_meanSC = [weekly_meanSC; mean(thisweek_SC(:,2))];
%             weekly_meanwear = [weekly_meanwear; mean(thisweek_SC(:,3))];
%         end
%     end
%     
%     this_weekly_raw(i).norm_SC(1).weekly_meanSC = weekly_meanSC;
%     this_weekly_raw(i).norm_SC(1).weekly_meanwear = weekly_meanwear;
%     
%     norm_factor = [];
%     for j = 1:length(weekly_meanwear)
%         norm_factor = [norm_factor; maxwear/weekly_meanwear(j)];
%     end
%     ydata = weekly_meanSC.*norm_factor;
%     this_weekly_raw(i).norm_SC(1).weekly_norm_meanSC = ydata;
%     
% end
% 
% weekly_raw_healthy = this_weekly_raw; 


%% Export data (optional)

this_weekly_raw = weekly_raw_healthy;

% Create table that stores significant variables for exporting.
varTypes = {'string','double','double','double','double'}; % Look at field names to determine var type.
fieldnames = {'Name','avg_end_expdist_mu','avg_int_sigma','avg_freq','product'};
T = table('Size',[length(this_weekly_raw) length(varTypes)],'VariableTypes',varTypes,'VariableNames',fieldnames);

names = {}; endurance = []; intensity = []; freq = []; product = [];
for i = 1:length(this_weekly_raw)
    names = [names; this_weekly_raw(i).name];
    intensity = [intensity; this_weekly_raw(i).avg_int_sigma];
    if isempty(this_weekly_raw(i).avg_end_expdist_mu)
        endurance = [endurance; 0];
        freq = [freq; 0];
    else
        endurance = [endurance; this_weekly_raw(i).avg_end_expdist_mu];
        freq = [freq; this_weekly_raw(i).avg_freq];
    end
    product = [product; this_weekly_raw(i).avg_end_expdist_mu*this_weekly_raw(i).avg_int_sigma*this_weekly_raw(i).avg_freq];
end

T(:,1) = names;
T(:,2) = num2cell(endurance);
T(:,3) = num2cell(intensity);
T(:,4) = num2cell(freq);
T(:,5) = num2cell(product);

delete 'sum_stats.xlsx'
writetable(T, 'sum_stats.xlsx')


%% Graph: endurance/frequency vs intensity/frequency for healthy subjects

% Run previous section of code before running this one.

%%% Figure 1
xdata = intensity./freq;
ydata = endurance./freq;
figure(1)
hold on
scatter(xdata,ydata)
title('Endurance/freq vs Intensity/freq for Healthy Subjects')
xlabel('Intensity / Frequency')
ylabel('Endurance / Frequency')

b = num2str(healthy_indices); c = cellstr(b);
dx = 0.01; dy = 0.05; % displacement so the text does not overlay the data points
text(xdata+dx, ydata, c);
hold off


%%% Figure 2
xdata = log10(intensity./freq);
ydata = log10(endurance./freq);
figure(2)
hold on
scatter(xdata,ydata)
[p1,S1] = polyfit(xdata, ydata, 1);
y = @(x) (x*p1(1) + p1(2));
fplot(y,'-k','LineWidth',1)
title('log(endurance/freq) vs log(intensity/freq) for Healthy Subjects')
xlabel('log(Intensity / Frequency)')
ylabel('log(Endurance / Frequency)')
xlim([-0.3 0.5])
ylim([-1.1 -0.4])
TextLocation(['y = ' num2str(sprintf('%.3g',p1(1))), '*x + ' num2str(sprintf('%.3g',p1(2))),...
    newline,newline,'Norm of residuals = ' num2str(sprintf('%.3g',S1.normr))],...
    'Location','northwest');

b = num2str(healthy_indices); c = cellstr(b);
dx = 0.02; dy = 0.02; % displacement so the text does not overlay the data points
text(xdata-dx, ydata+dy, c);
hold off


%%% Figure 3
figure(3)
hold on
scatter(intensity,endurance)
title('Endurance vs Intensity for Healthy Subjects')
xlabel('Intensity')
ylabel('Endurance')

b = num2str(healthy_indices); c = cellstr(b);
dx = 0.1; dy = 0.05; % displacement so the text does not overlay the data points
text(intensity+dx, endurance, c);
hold off


%% Find # of minutes of sustained walking for COPD patients
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

this_weekly_raw = weekly_raw_COPD;

% Gather # of minutes of sustained walking & avg SPM on week-by-week basis
for i = 1:length(this_weekly_raw)
    this_weekly_raw(i).sustained = struct([]);
    
    for j = 1:length(this_weekly_raw(i).weekly_raw)
        this_weekly_raw(i).sustained(j).WeekNo = this_weekly_raw(i).weekly_raw(j).WeekNo;
        this_weekly_raw(i).sustained(j).YearNo = this_weekly_raw(i).weekly_raw(j).YearNo;
        this_weekly_raw(i).sustained(j).duration = [];
        this_weekly_raw(i).sustained(j).avgSPM = [];
        this_weekly_raw(i).sustained(j).rawdata = struct([]);
        
        this_week = this_weekly_raw(i).weekly_raw(j);     
        this_time = this_week.all_time;
 
        time_start = datevec(datenum(this_time(1))); % Time at start of week
        time_end = datevec(datenum(this_time)); % Time at end of week
        elapsed_time = etime(time_end,time_start); % in seconds
 
        this_SC = this_week.all_SC;
        this_HR = this_week.all_HR; % Unused for now (2/25/2021).
        % [h,m,s] = hms(this_time);
        
        suswalk_idx = [];
        
        for k = 1:length(this_SC)-1
            if this_SC(k) > SCthresh
                streak = 1;
                while k+streak < length(this_SC)+1
                    if this_SC(k+streak) > SCthresh
                        if elapsed_time(k+streak) - elapsed_time(k) == streak*60
                            streak = streak + 1;
                        else
                            break
                        end
                    else
                        break
                    end
                end
                
                if k+streak > length(this_SC) % For edge case where sustained walking occurs at the very end of vector
                    streak = streak - 1;
                end
                
                if streak == 1
                    continue
                else
                    suswalk_idx = [suswalk_idx; k];
                    this_weekly_raw(i).sustained(j).duration = [this_weekly_raw(i).sustained(j).duration; streak];
                    this_weekly_raw(i).sustained(j).avgSPM = [this_weekly_raw(i).sustained(j).avgSPM; mean(this_SC(k:k+streak-1))];
                end
            end
        end
        
        % Remove "streaks" of sustained walking that are within the same ambulation.
        final_duration = [];
        final_avgSPM = [];
        final_suswalk_idx = [];
        idx = 1;
        while idx < length(this_weekly_raw(i).sustained(j).duration)+1
            this_duration = this_weekly_raw(i).sustained(j).duration(idx);
            final_duration = [final_duration; this_duration];
            final_avgSPM = [final_avgSPM; this_weekly_raw(i).sustained(j).avgSPM(idx+this_duration-2)];
            final_suswalk_idx = [final_suswalk_idx; suswalk_idx(idx)];
            idx = idx + this_weekly_raw(i).sustained(j).duration(idx) - 1;
        end
        this_weekly_raw(i).sustained(j).final_duration = final_duration; 
        this_weekly_raw(i).sustained(j).final_avgSPM = final_avgSPM;
        this_weekly_raw(i).sustained(j).compileddata = horzcat(final_duration, final_avgSPM);
        
        % Calculate 1/e value of weekly duration histogram
        if ~isempty(this_weekly_raw(i).sustained(j).compileddata)
            dur_data = this_weekly_raw(i).sustained(j).compileddata(:,1);
            dur_data_mod = [];
            for k = 1:length(dur_data)
                if dur_data(k) < min_dur
                else
                    dur_data_mod = [dur_data_mod; dur_data(k)];
                end
            end

            if length(unique(dur_data_mod)) > 1
                this_weekly_raw(i).sustained(j).dist_var = fitdist(dur_data_mod-min_dur,'exponential');
                this_weekly_raw(i).sustained(j).expdist_mu = this_weekly_raw(i).sustained(j).dist_var(1).mu+min_dur;
            end

            for k = 1:length(final_suswalk_idx)
                % Store raw data
                start_idx = final_suswalk_idx(k);
                end_idx = start_idx + final_duration(k) - 1;
                this_weekly_raw(i).sustained(j).rawdata(k).rawSPM = this_week.all_SC(start_idx:end_idx);
                this_weekly_raw(i).sustained(j).rawdata(k).rawHR = this_week.all_HR(start_idx:end_idx);
            end
        end
    end
    
    temp1 = []; temp2 = [];
    for j = 1:length(this_weekly_raw(i).sustained)
        if isfield(this_weekly_raw(i).sustained, 'expdist_mu')
            temp1 = [temp1; this_weekly_raw(i).sustained(j).expdist_mu];
            temp2 = [temp2; length(this_weekly_raw(i).sustained(j).final_duration)];
        end
    end
    this_weekly_raw(i).avg_end_expdist_mu = mean(temp1);
    this_weekly_raw(i).avg_freq = mean(temp2);
end
 
% Compile # of minutes of sustained walking & avg SPM across all weeks
for i = 1:length(this_weekly_raw)
    this_weekly_raw(i).compileddata = [];
    temp_duration = []; temp_avgSPM = [];
    
    for j = 1:length(this_weekly_raw(i).sustained)
        temp_duration = [temp_duration; this_weekly_raw(i).sustained(j).final_duration];
        temp_avgSPM = [temp_avgSPM; this_weekly_raw(i).sustained(j).final_avgSPM];
    end
    
    this_weekly_raw(i).compileddata = horzcat(temp_duration, temp_avgSPM);
    
end
 
% Find total number of sustained steps for each week
for i = 1:length(this_weekly_raw)
    for j = 1:length(this_weekly_raw(i).sustained)
        temp = this_weekly_raw(i).sustained(j).final_duration .* this_weekly_raw(i).sustained(j).final_avgSPM;
        this_weekly_raw(i).sustained(j).totalsussteps = sum(temp);
        this_weekly_raw(i).sustained(j).totalsteps = sum(this_weekly_raw(i).weekly_raw(j).all_SC);
        this_weekly_raw(i).sustained(j).susovertotal = this_weekly_raw(i).sustained(j).totalsussteps/this_weekly_raw(i).sustained(j).totalsteps;
    end
end
 
weekly_raw_COPD = this_weekly_raw;
 
 
%% Statistical Analysis of 3D Histogram/Heat Map

this_weekly_raw = weekly_raw_COPD;
 
for i = 1:length(this_weekly_raw)
    % this_weekly_raw(i).totalsteps_week = mean(this_weekly_raw(i).norm_SC(1).weekly_norm_meanSC)*7;
    
    temp1 = []; temp2 = []; temp3 = [];
    for j = 1:length(this_weekly_raw(i).sustained)
        temp1 = [temp1; this_weekly_raw(i).sustained(j).totalsussteps];
        temp2 = [temp2; this_weekly_raw(i).sustained(j).totalsteps];
        temp3 = [temp3; this_weekly_raw(i).sustained(j).susovertotal];
    end
    this_weekly_raw(i).avg_totalsussteps = mean(temp1);
    this_weekly_raw(i).avg_totalsteps = mean(temp2);
    this_weekly_raw(i).avg_susovertotal = mean(temp3);
end
 
weekly_raw_COPD = this_weekly_raw;
 
 
%% Graph heat map/3D histogram of sustained walking data for COPD subjects
 
this_weekly_raw = weekly_raw_COPD;

for i = 1:length(this_weekly_raw)
    if size(this_weekly_raw(i).compileddata,1) > histdatathresh
        % Specify centers of histogram bins
        bins_duration_max = max(this_weekly_raw(i).compileddata(:,1));
        % bins_duration_max = prctile(this_weekly_raw(i).compileddata(:,1),95);
        bins_avgSPM_max = max(this_weekly_raw(i).compileddata(:,2));
        bins_duration = 2:1:bins_duration_max;
        bins_duration_2D = 2:1:(bins_duration_max+1);
        bins_duration_2D = bins_duration_2D - 0.5;
        bins_avgSPM = 60:10:bins_avgSPM_max;
        
        [this_weekly_raw(i).N, this_weekly_raw(i).C] = hist3(this_weekly_raw(i).compileddata,'Ctrs',{bins_duration bins_avgSPM});
 
%         %%% Plot 3D histogram (sustained walking duration & avg SPM)
%         figure(i)
%         hold on
%         [this_weekly_raw(i).histcount,this_weekly_raw(i).histcenter] = hist3(this_weekly_raw(i).compileddata,'Ctrs',{bins_duration bins_avgSPM});
%         hist3(this_weekly_raw(i).compileddata,'Ctrs',{bins_duration bins_avgSPM},'CdataMode','auto')
%         title(sprintf('Histogram: Sustained Walking Duration & Avg SPM for %s',string(this_weekly_raw(i).name)))
%         xlabel('Duration (min)')
%         ylabel('SPM')
%         % colorbar
%         view(3)
%         TextLocation(['Avg # of sustained steps/week = ' num2str(sprintf('%d',round(this_weekly_raw(i).avg_totalsteps)))],...
%             'Location','best');
%         hold off
% 
%         %%% Plot heat map (sustained walking duration & avg SPM)
%         figure(i+length(this_weekly_raw))
%         hold on
%         hist3(this_weekly_raw(i).compileddata,'Ctrs',{bins_duration bins_avgSPM},'CdataMode','auto')
%         title(sprintf('Heat Map: Sustained Walking Duration & Avg SPM for %s',string(this_weekly_raw(i).name)))
%         xlabel('Duration (min)')
%         ylabel('SPM')
%         colorbar
%         hold off
%         
%         %%% Plot 2D duration histogram (disregard avg SPM)
%         figure(i+2*length(this_weekly_raw))
%         hold on
%         histogram(this_weekly_raw(i).compileddata(:,1),bins_duration_2D)
%         % Statistics
%         dur_median = median(this_weekly_raw(i).compileddata(:,1));
%         dur_sig = std(this_weekly_raw(i).compileddata(:,1));
%         dur_sk = skewness(this_weekly_raw(i).compileddata(:,1));
%         [~,~,dur_ks,~] = kstest(this_weekly_raw(i).compileddata(:,1),'CDF',fitdist(this_weekly_raw(i).compileddata(:,1),'lognormal'));
%         TextLocation(['Median = ' num2str(sprintf('%.3g',dur_median)), ', SD = ' num2str(sprintf('%.3g',dur_sig)),...
%             newline, 'Skewness = ' num2str(sprintf('%.3g',dur_sk)), ', KS-statistics = ' num2str(sprintf('%.3g',dur_ks))],...
%             'Location','northeast');
%         title(sprintf('Sustained Walking Duration for %s',string(this_weekly_raw(i).name)))
%         xlabel('Duration (min)')
%         ylabel('Frequency')   
%         hold off
        
        %%% Plot modified 2D duration histogram: require minimum # of sustained minutes
        dur_data = this_weekly_raw(i).compileddata(:,1);
        dur_data_mod = [];
        for j = 1:length(dur_data)
            if dur_data(j) < min_dur
            else
                dur_data_mod = [dur_data_mod; dur_data(j)];
            end
        end
        
        if length(unique(dur_data_mod)) > 1 % arbitrary cutoff
            % ~isempty(dur_data_mod)
            figure(i+3*length(this_weekly_raw))
            hold on
            histfit(dur_data_mod-min_dur,length(bins_duration_2D),'exponential')
            % histogram(dur_data_mod,length(bins_duration_2D))
            this_weekly_raw(i).end_dist_var = fitdist(dur_data_mod-min_dur,'exponential');
            
            % Statistics
            dur_median = median(dur_data_mod);
            dur_sig = std(dur_data_mod);
            dur_sk = skewness(dur_data_mod);
            if length(unique(dur_data_mod)) < 2
                dur_ks = NaN;
            else
                [~,~,dur_ks,~] = kstest(dur_data_mod,'CDF',fitdist(dur_data_mod,'lognormal'));
            end
            % 95th percentile of exponential distribution (PDF)
            ninetyfive = -this_weekly_raw(i).end_dist_var(1).mu * log(0.05);
            
            TextLocation(['Median = ' num2str(sprintf('%.3g',dur_median)), ', SD = ' num2str(sprintf('%.3g',dur_sig)),...
                newline, 'Skewness = ' num2str(sprintf('%.3g',dur_sk)), ', KS-statistics = ' num2str(sprintf('%.3g',dur_ks)),...
                newline, newline, '1/e value = ' num2str(sprintf('%.3g',this_weekly_raw(i).end_dist_var(1).mu)),...
                ', 95th percentile = ' num2str(sprintf('%.3g',ninetyfive))],...
                'Location','northeast');
            title(sprintf('Sustained Walking Duration for %s',string(this_weekly_raw(i).name)))
            xlabel('Duration (min)')
            ylabel('Frequency')   
            hold off
        end
        
%         %%% Plot 2D avg SPM histogram (disregard duration)
%         %%% THIS DIDN'T REALLY WORK THOUGH
%         figure(i+4*length(this_weekly_raw))
%         hold on
%         this_data = this_weekly_raw(i).compileddata(:,2);
%         histogram(this_data)
%         % Statistics
%         med = median(this_data);
%         sig = std(this_data);
%         sk = skewness(this_data);
%         [~,~,ks,~] = kstest(this_data,'CDF',fitdist(this_data,'exponential'));
%         TextLocation(['Median = ' num2str(sprintf('%.3g',med)), ', SD = ' num2str(sprintf('%.3g',sig)),...
%             newline, 'Skewness = ' num2str(sprintf('%.3g',sk)), ', KS-statistics = ' num2str(sprintf('%.3g',ks))],...
%             'Location','northeast');
%         title(sprintf('Avg SPM during Sustained Walking for %s',string(this_weekly_raw(i).name)))
%         xlabel('SPM')
%         ylabel('Frequency')   
%         hold off
 
    end
 
end
 
weekly_raw_COPD = this_weekly_raw;
 
 
%% Graph scatter plot of endurance vs ambulation frequency for COPD subjects
 
this_weekly_raw = weekly_raw_COPD;

for i = 1:length(this_weekly_raw)
    xdata = []; ydata = [];
    if isfield(this_weekly_raw(i).sustained, 'expdist_mu')
        for j = 1:length(this_weekly_raw(i).sustained)
            if ~isempty(this_weekly_raw(i).sustained(j).expdist_mu)
                xdata = [xdata; length(this_weekly_raw(i).sustained(j).final_duration)];
                ydata = [ydata; this_weekly_raw(i).sustained(j).expdist_mu];
            end
        end

        figure(i)
        hold on
        scatter(xdata,ydata)
        [p,S] = polyfit(xdata, ydata, 1);
        y = @(x) (x*p(1) + p(2));
        fplot(y,'-k','LineWidth',1.5)
        title(sprintf('Endurance (1/e) vs Ambulation Frequency for %s',string(this_weekly_raw(i).name)))
        xlabel('Ambulation Frequency')
        ylabel('Endurance (1/e of exp fit of duration histogram) (min)')
        TextLocation(['y = ' num2str(sprintf('%.3g',p(1))), '*x + ' num2str(sprintf('%.3g',p(2))),...
            newline,newline,'Norm of residuals = ' num2str(sprintf('%.3g',S.normr))],...
            'Location','northeast');
        hold off
    end
end
 
 
%% SC Histogram: "Intensity" Parameter for COPD patients (exponential fit)
 
% this_weekly_raw = weekly_raw_COPD;
%  
% for i = 1:length(this_weekly_raw)
%     this_weekly_raw(i).SChist = struct([]);
%     
%     for j = 1:length(this_weekly_raw(i).weekly_raw)
%         this_weekly_raw(i).SChist(j).WeekNo = this_weekly_raw(i).weekly_raw(j).WeekNo;
%         this_weekly_raw(i).SChist(j).YearNo = this_weekly_raw(i).weekly_raw(j).YearNo;
%         this_weekly_raw(i).SChist(j).dist_var = fitdist(this_weekly_raw(i).weekly_raw(j).SC,'exponential');
%         this_weekly_raw(i).SChist(j).expdist_mu = this_weekly_raw(i).SChist(j).dist_var(1).mu;
%     end
%     
%     temp1 = []; temp2 = [];
%     for j = 1:length(this_weekly_raw(i).SChist)
%         temp1 = [temp1; this_weekly_raw(i).SChist(j).expdist_mu];
%         temp2 = [temp2; this_weekly_raw(i).weekly_raw(j).SC];
%     end
%     this_weekly_raw(i).avg_int_expdist_mu = mean(temp1);
%     this_weekly_raw(i).std_int_expdist_mu = std(temp1);
%     
% %     figure(i)
% %     hold on
% %     histfit(temp2,max(temp2),'exponential')
% %     % Statistics
% %     med = median(temp2);
% %     sig = std(temp2);
% %     sk = skewness(temp2);
% %     if length(unique(temp2)) < 2
% %         ks = NaN;
% %     else
% %         [~,~,ks,~] = kstest(temp2,'CDF',fitdist(temp2,'exponential'));
% %     end
% %                
% %     TextLocation(['Median = ' num2str(sprintf('%.3g',med)), ', SD = ' num2str(sprintf('%.3g',sig)),...
% %         newline, 'Skewness = ' num2str(sprintf('%.3g',sk)), ', KS-statistics = ' num2str(sprintf('%.3g',ks)),...
% %         newline, newline, 'Avg 1/e value = ' num2str(sprintf('%.3g',this_weekly_raw(i).avg_expdist_mu))],...
% %         'Location','northeast');
% %     title(sprintf('Consolidated SC Histogram for %s',string(this_weekly_raw(i).name)))
% %     xlabel('SPM')
% %     ylabel('Frequency')   
% %     hold off
%     
% end
%  
% weekly_raw_COPD = this_weekly_raw;


%% SC Histogram: "Intensity" Parameter for COPD patients (lognormal fit)

this_weekly_raw = weekly_raw_COPD;

for i = 1:length(this_weekly_raw)
    this_weekly_raw(i).SChist = struct([]);
    
    for j = 1:length(this_weekly_raw(i).weekly_raw)
        this_weekly_raw(i).SChist(j).WeekNo = this_weekly_raw(i).weekly_raw(j).WeekNo;
        this_weekly_raw(i).SChist(j).YearNo = this_weekly_raw(i).weekly_raw(j).YearNo;
        this_weekly_raw(i).SChist(j).dist_var = fitdist(this_weekly_raw(i).weekly_raw(j).SC,'lognormal');
        this_weekly_raw(i).SChist(j).lndist_mu = this_weekly_raw(i).SChist(j).dist_var(1).mu;
        % this_weekly_raw(i).SChist(j).lndist_sigma = this_weekly_raw(i).SChist(j).dist_var(1).sigma;
    end
    
    temp1 = []; temp2 = []; temp3 = [];
    for j = 1:length(this_weekly_raw(i).SChist)
        temp1 = [temp1; this_weekly_raw(i).SChist(j).lndist_mu];
        temp2 = [temp2; this_weekly_raw(i).weekly_raw(j).SC];
        % temp3 = [temp3; this_weekly_raw(i).SChist(j).lndist_sigma];
    end
    this_weekly_raw(i).avg_int_lndist_mu = mean(temp1);
    this_weekly_raw(i).std_int_lndist_mu = std(temp1);
    % this_weekly_raw(i).avg_int_lndist_sigma = mean(temp3);
    
    figure(i)
    hold on
    histfit(temp2,max(temp2),'lognormal')
    pd = fitdist(temp2,'lognormal');
    % Statistics
    med = median(temp2);
    sig = std(temp2);
    this_weekly_raw(i).avg_int_sigma = sig; % Save this into structure
    sk = skewness(temp2);
    if length(unique(temp2)) < 2
        ks = NaN;
    else
        [~,~,ks,~] = kstest(temp2,'CDF',fitdist(temp2,'lognormal'));
    end
               
    TextLocation(['Median = ' num2str(sprintf('%.3g',med)), ', SD = ' num2str(sprintf('%.3g',sig)),...
        newline, 'Skewness = ' num2str(sprintf('%.3g',sk)), ', KS-statistics = ' num2str(sprintf('%.3g',ks))],...
        'Location','northeast');
    title(sprintf('Consolidated SC Histogram for %s',string(this_weekly_raw(i).name)))
    xlabel('SPM')
    ylabel('Frequency')   
    hold off
    
end

weekly_raw_COPD = this_weekly_raw;


%% Optional: Condense SC data for COPD subjects

this_weekly_raw = weekly_raw_COPD;
 
for i = 1:length(this_weekly_raw)
    temp1 = [];
    for j = 1:length(this_weekly_raw(i).weekly_raw)
        temp1 = [temp1; this_weekly_raw(i).weekly_raw(j).SC];
    end
    this_weekly_raw(i).allSC = temp1;
end

weekly_raw_COPD = this_weekly_raw;

figure(1)
histfit


%% Graph scatter plot of intensity vs ambulation frequency for COPD subjects
 
this_weekly_raw = weekly_raw_COPD;
 
for i = 1:length(this_weekly_raw)
    xdata = []; ydata = [];
    for j = 1:length(this_weekly_raw(i).SChist)
        if isfield(this_weekly_raw(i).sustained, 'expdist_mu')
            if ~isempty(this_weekly_raw(i).sustained(j).expdist_mu)
                xdata = [xdata; this_weekly_raw(i).sustained(j).expdist_mu];
                ydata = [ydata; this_weekly_raw(i).SChist(j).expdist_mu];
            end
        end
    end
    
    if ~isempty(xdata)
        figure(i)
        hold on
        scatter(xdata,ydata)
        [p,S] = polyfit(xdata, ydata, 1);
        y = @(x) (x*p(1) + p(2));
        fplot(y,'-k','LineWidth',1.5)
        title(sprintf('Intensity (1/e) vs Endurance (1/e) for %s',string(this_weekly_raw(i).name)))
        xlabel('Endurance (1/e of exp fit of SC histogram) (min)')
        ylabel('Intensity (1/e of exp fit of SC histogram) (SPM)')
        TextLocation(['y = ' num2str(sprintf('%.3g',p(1))), '*x + ' num2str(sprintf('%.3g',p(2))),...
            newline,newline,'Norm of residuals = ' num2str(sprintf('%.3g',S.normr))],...
            'Location','northeast');
        hold off
    end
end
 
 
%% Graph scatter plot of intensity vs ambulation frequency for COPD subjects
 
this_weekly_raw = weekly_raw_COPD;
 
for i = 1:length(this_weekly_raw)
    xdata = []; ydata = [];
    for j = 1:length(this_weekly_raw(i).SChist)
        if ~isempty(this_weekly_raw(i).SChist(j).expdist_mu)
            xdata = [xdata; length(this_weekly_raw(i).sustained(j).final_duration)];
            ydata = [ydata; this_weekly_raw(i).SChist(j).expdist_mu];
        end
    end
    figure(i)
    hold on
    scatter(xdata,ydata)
    [p,S] = polyfit(xdata, ydata, 1);
    y = @(x) (x*p(1) + p(2));
    fplot(y,'-k','LineWidth',1.5)
    title(sprintf('Intensity (1/e) vs Ambulation Frequency for %s',string(this_weekly_raw(i).name)))
    xlabel('Ambulation Frequency')
    ylabel('Intensity (1/e of exp fit of SC histogram) (SPM)')
    TextLocation(['y = ' num2str(sprintf('%.3g',p(1))), '*x + ' num2str(sprintf('%.3g',p(2))),...
        newline,newline,'Norm of residuals = ' num2str(sprintf('%.3g',S.normr))],...
        'Location','northeast');
    hold off
end
 
 
%% Find avg steps per day (using data between 6am - 10pm) for COPD Subjects
 
% this_weekly_raw = weekly_raw_COPD;
% 
% % If patient wore the Fitbit the entire time between 6am-10pm, how many data
% % points is that?
% maxwear = 60*16; % minutes
% 
% for i = 1:length(this_weekly_raw)
%     this_weekly_raw(i).norm_SC = struct([]);
%     week = []; year = [];
%     ydata = [];
%     
%     weekly_meanSC = [];
%     weekly_meanwear = [];
%     
%     for j = 1:length(this_weekly_raw(i).weekly_raw)
%         thisweek_days = [];
%         if isempty(this_weekly_raw(i).weekly_raw(j).all_time)
%         else
%             week = [week; this_weekly_raw(i).weekly_raw(j).WeekNo];
%             year = [year; this_weekly_raw(i).weekly_raw(j).YearNo];
%             for k = 1:length(this_weekly_raw(i).weekly_raw(j).all_time)
%                 thisweek_days = [thisweek_days; day(this_weekly_raw(i).weekly_raw(j).all_time(k))];
%             end
%             thisweek_uniquedays = unique(thisweek_days);
%             thisweek_SC = zeros(length(thisweek_uniquedays),3);
%             for k = 1:length(thisweek_uniquedays)
%                 fulldayentries = find(day(this_weekly_raw(i).weekly_raw(j).all_time) == ...
%                     thisweek_uniquedays(k));
%                 
%                 daytimeentries = [];
%                 daytimeSC = [];
%                 
%                 for l = fulldayentries(1):fulldayentries(length(fulldayentries))
%                     if hour(this_weekly_raw(i).weekly_raw(j).all_time(l)) > 5 &&... 
%                             hour(this_weekly_raw(i).weekly_raw(j).all_time(l)) < 22
%                         daytimeentries = [daytimeentries; l];
%                     else
%                     end
%                 end
%                 for l = 1:length(daytimeentries)
%                     daytimeSC = [daytimeSC; this_weekly_raw(i).weekly_raw(j).all_SC(daytimeentries(l))];
%                 end
%                 total_daytimeSC = sum(daytimeSC);
%                 total_daytime_weartime = length(daytimeSC);
%                 thisweek_SC(k,1) = thisweek_uniquedays(k);
%                 thisweek_SC(k,2) = total_daytimeSC;
%                 thisweek_SC(k,3) = total_daytime_weartime;
%             end
%             weekly_meanSC = [weekly_meanSC; mean(thisweek_SC(:,2))];
%             weekly_meanwear = [weekly_meanwear; mean(thisweek_SC(:,3))];
%         end
%     end
%     
%     this_weekly_raw(i).norm_SC(1).weekly_meanSC = weekly_meanSC;
%     this_weekly_raw(i).norm_SC(1).weekly_meanwear = weekly_meanwear;
%     
%     norm_factor = [];
%     for j = 1:length(weekly_meanwear)
%         norm_factor = [norm_factor; maxwear/weekly_meanwear(j)];
%     end
%     ydata = weekly_meanSC.*norm_factor;
%     this_weekly_raw(i).norm_SC(1).weekly_norm_meanSC = ydata;
%     
% end
% 
% weekly_raw_COPD = this_weekly_raw; 
 
 
%% Export data (optional)
 
this_weekly_raw = weekly_raw_COPD;
 
% Create table that stores significant variables for exporting.
varTypes = {'string','double','double','double','double'}; % Look at field names to determine var type.
fieldnames = {'Name','avg_end_expdist_mu','avg_int_sigma','avg_freq','product'};
T = table('Size',[length(this_weekly_raw) length(varTypes)],'VariableTypes',varTypes,'VariableNames',fieldnames);

names = {}; endurance = []; intensity = []; freq = []; product = [];
for i = 1:length(this_weekly_raw)
    names = [names; this_weekly_raw(i).name];
    intensity = [intensity; this_weekly_raw(i).avg_int_sigma];
    if isempty(this_weekly_raw(i).avg_end_expdist_mu)
        endurance = [endurance; 0];
        freq = [freq; 0];
    else
        endurance = [endurance; this_weekly_raw(i).avg_end_expdist_mu];
        freq = [freq; this_weekly_raw(i).avg_freq];
    end
    product = [product; this_weekly_raw(i).avg_end_expdist_mu*this_weekly_raw(i).avg_int_sigma*this_weekly_raw(i).avg_freq];
end

T(:,1) = names;
T(:,2) = num2cell(endurance);
T(:,3) = num2cell(intensity);
T(:,4) = num2cell(freq);
T(:,5) = num2cell(product);

delete 'sum_stats.xlsx'
writetable(T, 'sum_stats.xlsx')


%% Graph: endurance/frequency vs intensity/frequency for COPD patients

% Run previous section of code before running this one.

xdata = intensity./freq;
ydata = endurance./freq;
figure(1)
hold on
scatter(xdata,ydata)
title('Endurance/freq vs Intensity/freq for COPD Patients')
xlabel('Intensity / Frequency')
ylabel('Endurance / Frequency')

a = [2;3;4;5;6;7;8;9;10;12;13;14;15;16;17;18;19;20;21;22];
b = num2str(a); c = cellstr(b);
dx = 0.1; dy = 0.05; % displacement so the text does not overlay the data points
text(xdata-dx, ydata+dy, c);
hold off

figure(2)
hold on
scatter(intensity,endurance)
title('Endurance vs Intensity for COPD Patients')
xlabel('Intensity')
ylabel('Endurance')

b = num2str(a); c = cellstr(b);
dx = 0.15; dy = 0.05; % displacement so the text does not overlay the data points
text(intensity+dx, endurance, c);
hold off
 
