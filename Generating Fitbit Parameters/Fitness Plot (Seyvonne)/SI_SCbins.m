%% Clear workspace

clear;clc;close all;

%% Organize data for healthy, PC, COPD, and PAH patients

load('Daily_Healthy_PC_110220.mat'); % This structure contains data for healthy, PC, and COPD patients.
load('daily_PAH_101220.mat'); % This structure contains data for PAH patients.

% List indices for healthy, PC, and COPD patients in structure
healthy_indices = [2; 3; 7; 9; 13];
PC_indices = [4; 6];
COPD_indices = [10; 11; 14; 15; 16];

% List indices for PAH patients in structure
PAH_ex_indices = [6, 16, 18, 24, 25, 29]; % Patient indices to exclude
PAH_indices = 1:1:30; % Array with all PAH patients (even patients to exclude)

% List of patient indices to be included in analysis
for i = 1:length(PAH_ex_indices)
    PAH_indices(PAH_indices == PAH_ex_indices(i)) = [];
end

% Create new structures for each patient population.
weekly_raw_healthy = struct([]);
weekly_raw_PC = struct([]);
weekly_raw_COPD = struct([]);
weekly_raw_PAH = struct([]);

% for each individual in daily_Healthy structure (excludes PAH patients)
for n = 1:length(daily_Healthy)
    if ismember(n,healthy_indices)
        weekly_raw_healthy = [weekly_raw_healthy, daily_Healthy(n).weekly_raw];
    elseif ismember(n,PC_indices)
        weekly_raw_PC = [weekly_raw_PC, daily_Healthy(n).weekly_raw];
    elseif ismember(n,COPD_indices)
        weekly_raw_COPD = [weekly_raw_COPD, daily_Healthy(n).weekly_raw];
    else
        continue;
    end
end

% for each individual in daily_PAH structure
for n = 1:length(daily_PAH)
    if ismember(n,PAH_indices)
        weekly_raw_PAH = [weekly_raw_PAH, daily_PAH(n).weekly_raw];
    else
        continue;
    end
end

%%% General Parameters
% # of SC data points needed in order to count
screq = 10; % Can change this. Minimum is 1 data point.

% Create matrix to store slope & intercept for linear best fit line.
lin_polyfit_vals = zeros(4,2);


%% Plot SC & HR against time (this didn't really work)

% entry = 5; % change this
% hold on
% plot(datenum(weekly_raw_COPD(entry).all_time), weekly_raw_COPD(entry).all_HR)
% bar(datenum(weekly_raw_COPD(entry).all_time), weekly_raw_COPD(entry).all_SC)
% hold off


%% Plot HR vs SC (for SC > 0) for Healthy Subjects

% Plot representative scatter plots of HR vs SC (for SC > 0) for each healthy individual.
repscat_healthy = zeros(length(healthy_indices),4);
for i = 1:length(healthy_indices)
    repscat_healthy(i,1) = healthy_indices(i);
    repscat_healthy(i,2) = length(daily_Healthy(healthy_indices(i)).weekly_raw);
end

repscat_healthy(:,3) = cumsum(repscat_healthy(:,2),1);

% Find the week with the most data points for plotting.
for i = 1:length(healthy_indices)
    if i == 1
        max_HRSC = 0;
        for j = 1:repscat_healthy(i,3)
            if length(weekly_raw_healthy(j).SC) > max_HRSC
                max_HRSC = length(weekly_raw_healthy(j).SC);
                repscat_healthy(i,4) = j;
            end
        end
    else
        max_HRSC = 0;
        for j = (repscat_healthy(i-1,3)+1):repscat_healthy(i,3)
            if length(weekly_raw_healthy(j).SC) > max_HRSC
                max_HRSC = length(weekly_raw_healthy(j).SC);
                repscat_healthy(i,4) = j;
            end
        end
    end
end

for i = 1:length(healthy_indices)
    figure(i)
    entry = repscat_healthy(i,4);
    hold on
    scatter(weekly_raw_healthy(entry).SC, weekly_raw_healthy(entry).HR)
    p = polyfit(weekly_raw_healthy(entry).SC, weekly_raw_healthy(entry).HR, 1);
    x1 = linspace(0, max(weekly_raw_healthy(entry).SC));
    y1 = polyval(p, x1);
    plot(x1,y1)
    TextLocation(['y = ' num2str(p(1)), 'x + ' num2str(p(2))],'Location','best');
    title(sprintf('Healthy Subject %s',string(i)))
    xlabel('SC (SPM)')
    ylabel('HR (BPM)')
    hold off
end


%% Plot HR vs SC (for SC > 0) for PC Patients

% Plot representative scatter plots of HR vs SC (for SC > 0) for each PC patient.
repscat_PC = zeros(length(PC_indices),4);
for i = 1:length(PC_indices)
    repscat_PC(i,1) = PC_indices(i);
    repscat_PC(i,2) = length(daily_Healthy(PC_indices(i)).weekly_raw);
end

repscat_PC(:,3) = cumsum(repscat_PC(:,2),1);

% Find the week with the most data points for plotting.
for i = 1:length(PC_indices)
    if i == 1
        max_HRSC = 0;
        for j = 1:repscat_PC(i,3)
            if length(weekly_raw_PC(j).SC) > max_HRSC
                max_HRSC = length(weekly_raw_PC(j).SC);
                repscat_PC(i,4) = j;
            end
        end
    else
        max_HRSC = 0;
        for j = (repscat_PC(i-1,3)+1):repscat_PC(i,3)
            if length(weekly_raw_PC(j).SC) > max_HRSC
                max_HRSC = length(weekly_raw_PC(j).SC);
                repscat_PC(i,4) = j;
            end
        end
    end
end

for i = 1:length(PC_indices)
    figure(i+length(healthy_indices))
    entry = repscat_PC(i,4);
    hold on
    scatter(weekly_raw_PC(entry).SC, weekly_raw_PC(entry).HR)
    p = polyfit(weekly_raw_PC(entry).SC, weekly_raw_PC(entry).HR, 1);
    x1 = linspace(0, max(weekly_raw_PC(entry).SC));
    y1 = polyval(p, x1);
    plot(x1,y1)
    TextLocation(['y = ' num2str(p(1)), 'x + ' num2str(p(2))],'Location','best');
    title(sprintf('PC Patient %s',string(i)))
    xlabel('SC (SPM)')
    ylabel('HR (BPM)')
    hold off
end


%% Plot HR vs SC (for SC > 0) for COPD Patients

% Plot representative scatter plots of HR vs SC (for SC > 0) for each COPD patient.
repscat_COPD = zeros(length(COPD_indices),4);
for i = 1:length(COPD_indices)
    repscat_COPD(i,1) = COPD_indices(i);
    repscat_COPD(i,2) = length(daily_Healthy(COPD_indices(i)).weekly_raw);
end

repscat_COPD(:,3) = cumsum(repscat_COPD(:,2),1);

% Find the week with the most data points for plotting.
for i = 1:length(COPD_indices)
    if i == 1
        max_HRSC = 0;
        for j = 1:repscat_COPD(i,3)
            if length(weekly_raw_COPD(j).SC) > max_HRSC
                max_HRSC = length(weekly_raw_COPD(j).SC);
                repscat_COPD(i,4) = j;
            end
        end
    else
        max_HRSC = 0;
        for j = (repscat_COPD(i-1,3)+1):repscat_COPD(i,3)
            if length(weekly_raw_COPD(j).SC) > max_HRSC
                max_HRSC = length(weekly_raw_COPD(j).SC);
                repscat_COPD(i,4) = j;
            end
        end
    end
end

for i = 1:length(COPD_indices)
    figure(i+length(healthy_indices)+length(PC_indices))
    entry = repscat_COPD(i,4);
    hold on
    scatter(weekly_raw_COPD(entry).SC, weekly_raw_COPD(entry).HR)
    p = polyfit(weekly_raw_COPD(entry).SC, weekly_raw_COPD(entry).HR, 1);
    x1 = linspace(0, max(weekly_raw_COPD(entry).SC));
    y1 = polyval(p, x1);
    plot(x1,y1)
    TextLocation(['y = ' num2str(p(1)), 'x + ' num2str(p(2))],'Location','best');
    title(sprintf('COPD Patient %s',string(i+1)))
    xlabel('SC (SPM)')
    ylabel('HR (BPM)')
    hold off
end


%% Plot HR vs SC (for SC > 0) for PAH Patients

% Plot representative scatter plots of HR vs SC (for SC > 0) for each PAH patient.
repscat_PAH = zeros(length(PAH_indices),5);
for i = 1:length(PAH_indices)
    repscat_PAH(i,1) = PAH_indices(i);
    repscat_PAH(i,2) = length(daily_PAH(PAH_indices(i)).weekly_raw);
end

repscat_PAH(:,3) = cumsum(repscat_PAH(:,2),1);

% Find the week with the most data points for plotting.
for i = 1:length(PAH_indices)
    if i == 1
        max_HRSC = 0;
        for j = 1:repscat_PAH(i,3)
            if length(weekly_raw_PAH(j).SC) > max_HRSC
                max_HRSC = length(weekly_raw_PAH(j).SC);
                repscat_PAH(i,4) = j;
                repscat_PAH(i,5) = length(weekly_raw_PAH(j).SC);
            end
        end
    else
        max_HRSC = 0;
        for j = (repscat_PAH(i-1,3)+1):repscat_PAH(i,3)
            if length(weekly_raw_PAH(j).SC) > max_HRSC
                max_HRSC = length(weekly_raw_PAH(j).SC);
                repscat_PAH(i,4) = j;
                repscat_PAH(i,5) = length(weekly_raw_PAH(j).SC);
            end
        end
    end
end

for i = 15
    % 1:length(PAH_indices)
    figure(i)
    % entry = repscat_PAH(i,4);
    entry = 246; % Remove this later.
    hold on
    scatter(weekly_raw_PAH(entry).SC, weekly_raw_PAH(entry).HR)
    p = polyfit(weekly_raw_PAH(entry).SC, weekly_raw_PAH(entry).HR, 1);
    x1 = linspace(0, max(weekly_raw_PAH(entry).SC));
    y1 = polyval(p, x1);
    plot(x1,y1)
    TextLocation(['y = ' num2str(p(1)), 'x + ' num2str(p(2))],'Location','best');
    title(sprintf('PAH Patient %s',string(PAH_indices(i))))
    xlabel('SC (SPM)')
    ylabel('HR (BPM)')
    hold off
end


%% Determine average HR in different SC bins for Healthy Subjects

% Calculate average HR for different SC bins for each healthy individual.
scbins_healthy = struct([]);
this_scbins = scbins_healthy;
this_weekly_raw = weekly_raw_healthy;

for i = 1:length(healthy_indices)
    this_scbins(i).healthy_indices = healthy_indices(i);
    this_scbins(i).num_entries = repscat_healthy(i,2);
    this_scbins(i).cum_entries = repscat_healthy(i,3);
end
 
for i = 1:length(healthy_indices)
    if i == 1
        
        sc_lowsteps = []; hr_lowsteps = [];
        sc_medsteps = []; hr_medsteps = [];
        sc_pursteps = []; hr_pursteps = [];
        sc_slowwalk = []; hr_slowwalk = [];
        sc_medwalk = []; hr_medwalk = [];
        sc_briskwalk = []; hr_briskwalk = [];
        sc_fastamb = []; hr_fastamb = [];
        sc_run = []; hr_run = [];
        
        this_scbins(i).sc_lowsteps = []; this_scbins(i).hr_lowsteps = [];
        this_scbins(i).sc_medsteps = []; this_scbins(i).hr_medsteps = [];
        this_scbins(i).sc_pursteps = []; this_scbins(i).hr_pursteps = [];
        this_scbins(i).sc_slowwalk = []; this_scbins(i).hr_slowwalk = [];
        this_scbins(i).sc_medwalk = []; this_scbins(i).hr_medwalk = [];
        this_scbins(i).sc_briskwalk = []; this_scbins(i).hr_briskwalk = [];
        this_scbins(i).sc_fastamb = []; this_scbins(i).hr_fastamb = [];
        this_scbins(i).sc_run = []; this_scbins(i).hr_run = [];
            
        for j = 1:this_scbins(i).cum_entries
            lowsteps = find(this_weekly_raw(j).SC > 0 & this_weekly_raw(j).SC < 20);
            medsteps = find(this_weekly_raw(j).SC > 19 & this_weekly_raw(j).SC < 40);
            pursteps = find(this_weekly_raw(j).SC > 39 & this_weekly_raw(j).SC < 60);
            slowwalk = find(this_weekly_raw(j).SC > 59 & this_weekly_raw(j).SC < 80);
            medwalk = find(this_weekly_raw(j).SC > 79 & this_weekly_raw(j).SC < 100);
            briskwalk = find(this_weekly_raw(j).SC > 99 & this_weekly_raw(j).SC < 120);
            fastamb = find(this_weekly_raw(j).SC > 119 & this_weekly_raw(j).SC < 150);
            run = find(this_weekly_raw(j).SC > 149 & this_weekly_raw(j).SC < 250);
            
            if length(lowsteps) > screq - 1
                for k = 1:length(lowsteps)
                    entry = lowsteps(k);
                    sc_lowsteps = [sc_lowsteps; this_weekly_raw(j).SC(entry)];
                    hr_lowsteps = [hr_lowsteps; this_weekly_raw(j).HR(entry)];
                end
                this_scbins(i).sc_lowsteps = [this_scbins(i).sc_lowsteps; sc_lowsteps];
                this_scbins(i).hr_lowsteps = [this_scbins(i).hr_lowsteps; hr_lowsteps];
            else
                this_scbins(i).meansc_lowsteps = NaN;
                this_scbins(i).meanhr_lowsteps = NaN;
            end
            
            if length(medsteps) > screq - 1
                for k = 1:length(medsteps)
                    entry = medsteps(k);
                    sc_medsteps = [sc_medsteps; this_weekly_raw(j).SC(entry)];
                    hr_medsteps = [hr_medsteps; this_weekly_raw(j).HR(entry)];
                end
                this_scbins(i).sc_medsteps = [this_scbins(i).sc_medsteps; sc_medsteps];
                this_scbins(i).hr_medsteps = [this_scbins(i).hr_medsteps; hr_medsteps];
            else
                this_scbins(i).meansc_medsteps = NaN;
                this_scbins(i).meanhr_medsteps = NaN;
            end
            
            if length(pursteps) > screq - 1
                for k = 1:length(pursteps)
                    entry = pursteps(k);
                    sc_pursteps = [sc_pursteps; this_weekly_raw(j).SC(entry)];
                    hr_pursteps = [hr_pursteps; this_weekly_raw(j).HR(entry)];
                end
                this_scbins(i).sc_pursteps = [this_scbins(i).sc_pursteps; sc_pursteps];
                this_scbins(i).hr_pursteps = [this_scbins(i).hr_pursteps; hr_pursteps];
            else
                this_scbins(i).meansc_pursteps = NaN;
                this_scbins(i).meanhr_pursteps = NaN;       
            end
            
            if length(slowwalk) > screq - 1
                for k = 1:length(slowwalk)
                    entry = slowwalk(k);
                    sc_slowwalk = [sc_slowwalk; this_weekly_raw(j).SC(entry)];
                    hr_slowwalk = [hr_slowwalk; this_weekly_raw(j).HR(entry)];
                end
                this_scbins(i).sc_slowwalk = [this_scbins(i).sc_slowwalk; sc_slowwalk];
                this_scbins(i).hr_slowwalk = [this_scbins(i).hr_slowwalk; hr_slowwalk];
            else
                this_scbins(i).meansc_slowwalk = NaN;
                this_scbins(i).meanhr_slowwalk = NaN;
            end
            
            if length(medwalk) > screq - 1
                for k = 1:length(medwalk)
                    entry = medwalk(k);
                    sc_medwalk = [sc_medwalk; this_weekly_raw(j).SC(entry)];
                    hr_medwalk = [hr_medwalk; this_weekly_raw(j).HR(entry)];
                end
                this_scbins(i).sc_medwalk = [this_scbins(i).sc_medwalk; sc_medwalk];
                this_scbins(i).hr_medwalk = [this_scbins(i).hr_medwalk; hr_medwalk];
            else
                this_scbins(i).meansc_medwalk = NaN;
                this_scbins(i).meanhr_medwalk = NaN;
            end
            
            if length(briskwalk) > screq - 1
                for k = 1:length(briskwalk)
                    entry = briskwalk(k);
                    sc_briskwalk = [sc_briskwalk; this_weekly_raw(j).SC(entry)];
                    hr_briskwalk = [hr_briskwalk; this_weekly_raw(j).HR(entry)];
                end
                this_scbins(i).sc_briskwalk = [this_scbins(i).sc_briskwalk; sc_briskwalk];
                this_scbins(i).hr_briskwalk = [this_scbins(i).hr_briskwalk; hr_briskwalk];
            else
                this_scbins(i).meansc_briskwalk = NaN;
                this_scbins(i).meanhr_briskwalk = NaN;
            end
            
            if length(fastamb) > screq - 1
                for k = 1:length(fastamb)
                    entry = fastamb(k);
                    sc_fastamb = [sc_fastamb; this_weekly_raw(j).SC(entry)];
                    hr_fastamb = [hr_fastamb; this_weekly_raw(j).HR(entry)];
                end
                this_scbins(i).sc_fastamb = [this_scbins(i).sc_fastamb; sc_fastamb];
                this_scbins(i).hr_fastamb = [this_scbins(i).hr_fastamb; hr_fastamb];
            else
                this_scbins(i).meansc_fastamb = NaN;
                this_scbins(i).meanhr_fastamb = NaN;
            end
            
            if length(run) > screq - 1
                for k = 1:length(run)
                    entry = run(k);
                    sc_run = [sc_run; this_weekly_raw(j).SC(entry)];
                    hr_run = [hr_run; this_weekly_raw(j).HR(entry)];
                end
                this_scbins(i).sc_run = [this_scbins(i).sc_run; sc_run];
                this_scbins(i).hr_run = [this_scbins(i).hr_run; hr_run];
            else
                this_scbins(i).meansc_run = NaN;
                this_scbins(i).meanhr_run = NaN;
            end
            
        end
        
        this_scbins(i).meansc_lowsteps = mean(this_scbins(i).sc_lowsteps);
        this_scbins(i).meanhr_lowsteps = mean(this_scbins(i).hr_lowsteps);
        this_scbins(i).meansc_medsteps = mean(this_scbins(i).sc_medsteps);
        this_scbins(i).meanhr_medsteps = mean(this_scbins(i).hr_medsteps);
        this_scbins(i).meansc_pursteps = mean(this_scbins(i).sc_pursteps);
        this_scbins(i).meanhr_pursteps = mean(this_scbins(i).hr_pursteps);
        this_scbins(i).meansc_slowwalk = mean(this_scbins(i).sc_slowwalk);
        this_scbins(i).meanhr_slowwalk = mean(this_scbins(i).hr_slowwalk);
        this_scbins(i).meansc_medwalk = mean(this_scbins(i).sc_medwalk);
        this_scbins(i).meanhr_medwalk = mean(this_scbins(i).hr_medwalk);
        this_scbins(i).meansc_briskwalk = mean(this_scbins(i).sc_briskwalk);
        this_scbins(i).meanhr_briskwalk = mean(this_scbins(i).hr_briskwalk);
        this_scbins(i).meansc_fastamb = mean(this_scbins(i).sc_fastamb);
        this_scbins(i).meanhr_fastamb = mean(this_scbins(i).hr_fastamb);
        this_scbins(i).meansc_run = mean(this_scbins(i).sc_run);
        this_scbins(i).meanhr_run = mean(this_scbins(i).hr_run);
        
    else
        
        sc_lowsteps = []; hr_lowsteps = [];
        sc_medsteps = []; hr_medsteps = [];
        sc_pursteps = []; hr_pursteps = [];
        sc_slowwalk = []; hr_slowwalk = [];
        sc_medwalk = []; hr_medwalk = [];
        sc_briskwalk = []; hr_briskwalk = [];
        sc_fastamb = []; hr_fastamb = [];
        sc_run = []; hr_run = [];
        
        this_scbins(i).sc_lowsteps = []; this_scbins(i).hr_lowsteps = [];
        this_scbins(i).sc_medsteps = []; this_scbins(i).hr_medsteps = [];
        this_scbins(i).sc_pursteps = []; this_scbins(i).hr_pursteps = [];
        this_scbins(i).sc_slowwalk = []; this_scbins(i).hr_slowwalk = [];
        this_scbins(i).sc_medwalk = []; this_scbins(i).hr_medwalk = [];
        this_scbins(i).sc_briskwalk = []; this_scbins(i).hr_briskwalk = [];
        this_scbins(i).sc_fastamb = []; this_scbins(i).hr_fastamb = [];
        this_scbins(i).sc_run = []; this_scbins(i).hr_run = [];
            
        for j = ((this_scbins(i-1).cum_entries)+1):this_scbins(i).cum_entries
            lowsteps = find(this_weekly_raw(j).SC > 0 & this_weekly_raw(j).SC < 20);
            medsteps = find(this_weekly_raw(j).SC > 19 & this_weekly_raw(j).SC < 40);
            pursteps = find(this_weekly_raw(j).SC > 39 & this_weekly_raw(j).SC < 60);
            slowwalk = find(this_weekly_raw(j).SC > 59 & this_weekly_raw(j).SC < 80);
            medwalk = find(this_weekly_raw(j).SC > 79 & this_weekly_raw(j).SC < 100);
            briskwalk = find(this_weekly_raw(j).SC > 99 & this_weekly_raw(j).SC < 120);
            fastamb = find(this_weekly_raw(j).SC > 119 & this_weekly_raw(j).SC < 150);
            run = find(this_weekly_raw(j).SC > 149 & this_weekly_raw(j).SC < 250);
            
            if length(lowsteps) > screq - 1
                for k = 1:length(lowsteps)
                    entry = lowsteps(k);
                    sc_lowsteps = [sc_lowsteps; this_weekly_raw(j).SC(entry)];
                    hr_lowsteps = [hr_lowsteps; this_weekly_raw(j).HR(entry)];
                end
                this_scbins(i).sc_lowsteps = [this_scbins(i).sc_lowsteps; sc_lowsteps];
                this_scbins(i).hr_lowsteps = [this_scbins(i).hr_lowsteps; hr_lowsteps];
            else
                this_scbins(i).meansc_lowsteps = NaN;
                this_scbins(i).meanhr_lowsteps = NaN;
            end
            
            if length(medsteps) > screq - 1
                for k = 1:length(medsteps)
                    entry = medsteps(k);
                    sc_medsteps = [sc_medsteps; this_weekly_raw(j).SC(entry)];
                    hr_medsteps = [hr_medsteps; this_weekly_raw(j).HR(entry)];
                end
                this_scbins(i).sc_medsteps = [this_scbins(i).sc_medsteps; sc_medsteps];
                this_scbins(i).hr_medsteps = [this_scbins(i).hr_medsteps; hr_medsteps];
            else
                this_scbins(i).meansc_medsteps = NaN;
                this_scbins(i).meanhr_medsteps = NaN;
            end
            
            if length(pursteps) > screq - 1
                for k = 1:length(pursteps)
                    entry = pursteps(k);
                    sc_pursteps = [sc_pursteps; this_weekly_raw(j).SC(entry)];
                    hr_pursteps = [hr_pursteps; this_weekly_raw(j).HR(entry)];
                end
                this_scbins(i).sc_pursteps = [this_scbins(i).sc_pursteps; sc_pursteps];
                this_scbins(i).hr_pursteps = [this_scbins(i).hr_pursteps; hr_pursteps];
            else
                this_scbins(i).meansc_pursteps = NaN;
                this_scbins(i).meanhr_pursteps = NaN;       
            end
            
            if length(slowwalk) > screq - 1
                for k = 1:length(slowwalk)
                    entry = slowwalk(k);
                    sc_slowwalk = [sc_slowwalk; this_weekly_raw(j).SC(entry)];
                    hr_slowwalk = [hr_slowwalk; this_weekly_raw(j).HR(entry)];
                end
                this_scbins(i).sc_slowwalk = [this_scbins(i).sc_slowwalk; sc_slowwalk];
                this_scbins(i).hr_slowwalk = [this_scbins(i).hr_slowwalk; hr_slowwalk];
            else
                this_scbins(i).meansc_slowwalk = NaN;
                this_scbins(i).meanhr_slowwalk = NaN;
            end
            
            if length(medwalk) > screq - 1
                for k = 1:length(medwalk)
                    entry = medwalk(k);
                    sc_medwalk = [sc_medwalk; this_weekly_raw(j).SC(entry)];
                    hr_medwalk = [hr_medwalk; this_weekly_raw(j).HR(entry)];
                end
                this_scbins(i).sc_medwalk = [this_scbins(i).sc_medwalk; sc_medwalk];
                this_scbins(i).hr_medwalk = [this_scbins(i).hr_medwalk; hr_medwalk];
            else
                this_scbins(i).meansc_medwalk = NaN;
                this_scbins(i).meanhr_medwalk = NaN;
            end
            
            if length(briskwalk) > screq - 1
                for k = 1:length(briskwalk)
                    entry = briskwalk(k);
                    sc_briskwalk = [sc_briskwalk; this_weekly_raw(j).SC(entry)];
                    hr_briskwalk = [hr_briskwalk; this_weekly_raw(j).HR(entry)];
                end
                this_scbins(i).sc_briskwalk = [this_scbins(i).sc_briskwalk; sc_briskwalk];
                this_scbins(i).hr_briskwalk = [this_scbins(i).hr_briskwalk; hr_briskwalk];
            else
                this_scbins(i).meansc_briskwalk = NaN;
                this_scbins(i).meanhr_briskwalk = NaN;
            end
            
            if length(fastamb) > screq - 1
                for k = 1:length(fastamb)
                    entry = fastamb(k);
                    sc_fastamb = [sc_fastamb; this_weekly_raw(j).SC(entry)];
                    hr_fastamb = [hr_fastamb; this_weekly_raw(j).HR(entry)];
                end
                this_scbins(i).sc_fastamb = [this_scbins(i).sc_fastamb; sc_fastamb];
                this_scbins(i).hr_fastamb = [this_scbins(i).hr_fastamb; hr_fastamb];
            else
                this_scbins(i).meansc_fastamb = NaN;
                this_scbins(i).meanhr_fastamb = NaN;
            end
            
            if length(run) > screq - 1
                for k = 1:length(run)
                    entry = run(k);
                    sc_run = [sc_run; this_weekly_raw(j).SC(entry)];
                    hr_run = [hr_run; this_weekly_raw(j).HR(entry)];
                end
                this_scbins(i).sc_run = [this_scbins(i).sc_run; sc_run];
                this_scbins(i).hr_run = [this_scbins(i).hr_run; hr_run];
            else
                this_scbins(i).meansc_run = NaN;
                this_scbins(i).meanhr_run = NaN;
            end
            
        end
        
        this_scbins(i).meansc_lowsteps = mean(this_scbins(i).sc_lowsteps);
        this_scbins(i).meanhr_lowsteps = mean(this_scbins(i).hr_lowsteps);
        this_scbins(i).meansc_medsteps = mean(this_scbins(i).sc_medsteps);
        this_scbins(i).meanhr_medsteps = mean(this_scbins(i).hr_medsteps);
        this_scbins(i).meansc_pursteps = mean(this_scbins(i).sc_pursteps);
        this_scbins(i).meanhr_pursteps = mean(this_scbins(i).hr_pursteps);
        this_scbins(i).meansc_slowwalk = mean(this_scbins(i).sc_slowwalk);
        this_scbins(i).meanhr_slowwalk = mean(this_scbins(i).hr_slowwalk);
        this_scbins(i).meansc_medwalk = mean(this_scbins(i).sc_medwalk);
        this_scbins(i).meanhr_medwalk = mean(this_scbins(i).hr_medwalk);
        this_scbins(i).meansc_briskwalk = mean(this_scbins(i).sc_briskwalk);
        this_scbins(i).meanhr_briskwalk = mean(this_scbins(i).hr_briskwalk);
        this_scbins(i).meansc_fastamb = mean(this_scbins(i).sc_fastamb);
        this_scbins(i).meanhr_fastamb = mean(this_scbins(i).hr_fastamb);
        this_scbins(i).meansc_run = mean(this_scbins(i).sc_run);
        this_scbins(i).meanhr_run = mean(this_scbins(i).hr_run);
            
    end
end
    
%     this_scbins(i).sdhr_lowsteps = std(hr_lowsteps);
%     this_scbins(i).sdhr_medsteps = std(hr_medsteps);
%     this_scbins(i).sdhr_pursteps = std(hr_pursteps);
%     this_scbins(i).sdhr_slowwalk = std(hr_slowwalk);
%     this_scbins(i).sdhr_medwalk = std(hr_medwalk);
%     this_scbins(i).sdhr_briskwalk = std(hr_briskwalk);
%     this_scbins(i).sdhr_fastamb = std(hr_fastamb);
%     this_scbins(i).sdhr_run = std(hr_run);

scbins_healthy = this_scbins;


%% Graph average HR in different SC bins for Healthy Subjects

allhealthy_scbins = zeros(length(healthy_indices),8); % Change based on # bins.
this_scbins = scbins_healthy;
 
figure(1)
allxdata = [];
allydata = [];
hold on
for i = 1:length(healthy_indices)
    xdata = [this_scbins(i).meansc_lowsteps, this_scbins(i).meansc_medsteps,...
        this_scbins(i).meansc_pursteps, this_scbins(i).meansc_slowwalk, ...
        this_scbins(i).meansc_medwalk, this_scbins(i).meansc_briskwalk, ...
        this_scbins(i).meansc_fastamb, this_scbins(i).meansc_run];
    ydata = [this_scbins(i).meanhr_lowsteps, this_scbins(i).meanhr_medsteps,...
        this_scbins(i).meanhr_pursteps, this_scbins(i).meanhr_slowwalk, ...
        this_scbins(i).meanhr_medwalk, this_scbins(i).meanhr_briskwalk, ...
        this_scbins(i).meanhr_fastamb, this_scbins(i).meanhr_run];
    allhealthy_scbins(i,:) = ydata;
    for j = 1:length(xdata)
        nanind = find(isnan(xdata));
        if isempty(nanind)
            xdata = xdata;
        else
            for k = 1:length(nanind)
                xdata(nanind(length(nanind)-k+1)) = [];
                ydata(nanind(length(nanind)-k+1)) = [];
            end
        end
    end
    allxdata = [allxdata, xdata];
    allydata = [allydata, ydata];
    plot(xdata,ydata,'-o','LineWidth',2)
%     err_std = [this_scbins(i).sdhr_lowsteps, this_scbins(i).sdhr_medwalk, ...
%         this_scbins(i).sdhr_pursteps, this_scbins(i).sdhr_slowwalk, ...
%         this_scbins(i).sdhr_medwalk, this_scbins(i).sdhr_briskwalk, ...
%         this_scbins(i).sdhr_fastamb, this_scbins(i).sdhr_run];
%     errorbar(xdata, ydata, err_std, '-o')
    p = polyfit(xdata,ydata, 1);
    this_scbins(i).slope = p(1);
    this_scbins(i).intercept = p(2);
end

% % Old version: using polyfit on ALL data (no differentiation b/t subjects)
% p = polyfit(allxdata,allydata, 1);
% lin_polyfit_vals(1,1) = p(1);
% lin_polyfit_vals(1,2) = p(2);
% x1 = linspace(0,max(allxdata));
% y1 = polyval(p, x1);
% plot(x1,y1,'Color','k','LineWidth',2)
% TextLocation(['y = ' num2str(p(1)), 'x + ' num2str(p(2))],'Location','southeast');

% New version: taking average slope & intercept to plot LSF line
slope = [];
int = [];
for i = 1:length(this_scbins)
    slope = [slope this_scbins(i).slope];
    int = [int this_scbins(i).intercept];
end
lin_polyfit_vals(1,1) = mean(slope);
lin_polyfit_vals(1,2) = mean(int);
x2 = linspace(0,max(allxdata));
y2 = @(x2) (mean(slope)*x2 + mean(int));
fplot(y2,'-k','LineWidth',2)
TextLocation(['y = ' num2str(mean(slope)),... 
    'x + ' num2str(mean(int))],'Location','southeast');

title('Mean HR at Various SC Thresholds for Healthy Subjects')
xlabel('SPM')
ylabel('Mean HR (BPM)')
legend('Incheol','Nicole','PeterX','PeterS','Seyvonne','Location','best')
hold off

scbins_healthy = this_scbins;

%% Determine average HR in different SC bins for PC Patients

% Calculate average HR for different SC bins for each PC patient.
scbins_PC = struct([]);
this_scbins = scbins_PC;
this_weekly_raw = weekly_raw_PC;

for i = 1:length(PC_indices)
    this_scbins(i).PC_indices = PC_indices(i);
    this_scbins(i).num_entries = repscat_PC(i,2);
    this_scbins(i).cum_entries = repscat_PC(i,3);
end
 
for i = 1:length(PC_indices)
    if i == 1
        
        sc_lowsteps = []; hr_lowsteps = [];
        sc_medsteps = []; hr_medsteps = [];
        sc_pursteps = []; hr_pursteps = [];
        sc_slowwalk = []; hr_slowwalk = [];
        sc_medwalk = []; hr_medwalk = [];
        sc_briskwalk = []; hr_briskwalk = [];
        sc_fastamb = []; hr_fastamb = [];
        sc_run = []; hr_run = [];
        
        this_scbins(i).sc_lowsteps = []; this_scbins(i).hr_lowsteps = [];
        this_scbins(i).sc_medsteps = []; this_scbins(i).hr_medsteps = [];
        this_scbins(i).sc_pursteps = []; this_scbins(i).hr_pursteps = [];
        this_scbins(i).sc_slowwalk = []; this_scbins(i).hr_slowwalk = [];
        this_scbins(i).sc_medwalk = []; this_scbins(i).hr_medwalk = [];
        this_scbins(i).sc_briskwalk = []; this_scbins(i).hr_briskwalk = [];
        this_scbins(i).sc_fastamb = []; this_scbins(i).hr_fastamb = [];
        this_scbins(i).sc_run = []; this_scbins(i).hr_run = [];
            
        for j = 1:this_scbins(i).cum_entries
            lowsteps = find(this_weekly_raw(j).SC > 0 & this_weekly_raw(j).SC < 20);
            medsteps = find(this_weekly_raw(j).SC > 19 & this_weekly_raw(j).SC < 40);
            pursteps = find(this_weekly_raw(j).SC > 39 & this_weekly_raw(j).SC < 60);
            slowwalk = find(this_weekly_raw(j).SC > 59 & this_weekly_raw(j).SC < 80);
            medwalk = find(this_weekly_raw(j).SC > 79 & this_weekly_raw(j).SC < 100);
            briskwalk = find(this_weekly_raw(j).SC > 99 & this_weekly_raw(j).SC < 120);
            fastamb = find(this_weekly_raw(j).SC > 119 & this_weekly_raw(j).SC < 150);
            run = find(this_weekly_raw(j).SC > 149 & this_weekly_raw(j).SC < 250);
            
            if length(lowsteps) > screq - 1
                for k = 1:length(lowsteps)
                    entry = lowsteps(k);
                    sc_lowsteps = [sc_lowsteps; this_weekly_raw(j).SC(entry)];
                    hr_lowsteps = [hr_lowsteps; this_weekly_raw(j).HR(entry)];
                end
                this_scbins(i).sc_lowsteps = [this_scbins(i).sc_lowsteps; sc_lowsteps];
                this_scbins(i).hr_lowsteps = [this_scbins(i).hr_lowsteps; hr_lowsteps];
            else
                this_scbins(i).meansc_lowsteps = NaN;
                this_scbins(i).meanhr_lowsteps = NaN;
            end
            
            if length(medsteps) > screq - 1
                for k = 1:length(medsteps)
                    entry = medsteps(k);
                    sc_medsteps = [sc_medsteps; this_weekly_raw(j).SC(entry)];
                    hr_medsteps = [hr_medsteps; this_weekly_raw(j).HR(entry)];
                end
                this_scbins(i).sc_medsteps = [this_scbins(i).sc_medsteps; sc_medsteps];
                this_scbins(i).hr_medsteps = [this_scbins(i).hr_medsteps; hr_medsteps];
            else
                this_scbins(i).meansc_medsteps = NaN;
                this_scbins(i).meanhr_medsteps = NaN;
            end
            
            if length(pursteps) > screq - 1
                for k = 1:length(pursteps)
                    entry = pursteps(k);
                    sc_pursteps = [sc_pursteps; this_weekly_raw(j).SC(entry)];
                    hr_pursteps = [hr_pursteps; this_weekly_raw(j).HR(entry)];
                end
                this_scbins(i).sc_pursteps = [this_scbins(i).sc_pursteps; sc_pursteps];
                this_scbins(i).hr_pursteps = [this_scbins(i).hr_pursteps; hr_pursteps];
            else
                this_scbins(i).meansc_pursteps = NaN;
                this_scbins(i).meanhr_pursteps = NaN;       
            end
            
            if length(slowwalk) > screq - 1
                for k = 1:length(slowwalk)
                    entry = slowwalk(k);
                    sc_slowwalk = [sc_slowwalk; this_weekly_raw(j).SC(entry)];
                    hr_slowwalk = [hr_slowwalk; this_weekly_raw(j).HR(entry)];
                end
                this_scbins(i).sc_slowwalk = [this_scbins(i).sc_slowwalk; sc_slowwalk];
                this_scbins(i).hr_slowwalk = [this_scbins(i).hr_slowwalk; hr_slowwalk];
            else
                this_scbins(i).meansc_slowwalk = NaN;
                this_scbins(i).meanhr_slowwalk = NaN;
            end
            
            if length(medwalk) > screq - 1
                for k = 1:length(medwalk)
                    entry = medwalk(k);
                    sc_medwalk = [sc_medwalk; this_weekly_raw(j).SC(entry)];
                    hr_medwalk = [hr_medwalk; this_weekly_raw(j).HR(entry)];
                end
                this_scbins(i).sc_medwalk = [this_scbins(i).sc_medwalk; sc_medwalk];
                this_scbins(i).hr_medwalk = [this_scbins(i).hr_medwalk; hr_medwalk];
            else
                this_scbins(i).meansc_medwalk = NaN;
                this_scbins(i).meanhr_medwalk = NaN;
            end
            
            if length(briskwalk) > screq - 1
                for k = 1:length(briskwalk)
                    entry = briskwalk(k);
                    sc_briskwalk = [sc_briskwalk; this_weekly_raw(j).SC(entry)];
                    hr_briskwalk = [hr_briskwalk; this_weekly_raw(j).HR(entry)];
                end
                this_scbins(i).sc_briskwalk = [this_scbins(i).sc_briskwalk; sc_briskwalk];
                this_scbins(i).hr_briskwalk = [this_scbins(i).hr_briskwalk; hr_briskwalk];
            else
                this_scbins(i).meansc_briskwalk = NaN;
                this_scbins(i).meanhr_briskwalk = NaN;
            end
            
            if length(fastamb) > screq - 1
                for k = 1:length(fastamb)
                    entry = fastamb(k);
                    sc_fastamb = [sc_fastamb; this_weekly_raw(j).SC(entry)];
                    hr_fastamb = [hr_fastamb; this_weekly_raw(j).HR(entry)];
                end
                this_scbins(i).sc_fastamb = [this_scbins(i).sc_fastamb; sc_fastamb];
                this_scbins(i).hr_fastamb = [this_scbins(i).hr_fastamb; hr_fastamb];
            else
                this_scbins(i).meansc_fastamb = NaN;
                this_scbins(i).meanhr_fastamb = NaN;
            end
            
            if length(run) > screq - 1
                for k = 1:length(run)
                    entry = run(k);
                    sc_run = [sc_run; this_weekly_raw(j).SC(entry)];
                    hr_run = [hr_run; this_weekly_raw(j).HR(entry)];
                end
                this_scbins(i).sc_run = [this_scbins(i).sc_run; sc_run];
                this_scbins(i).hr_run = [this_scbins(i).hr_run; hr_run];
            else
                this_scbins(i).meansc_run = NaN;
                this_scbins(i).meanhr_run = NaN;
            end
            
        end
        
        this_scbins(i).meansc_lowsteps = mean(this_scbins(i).sc_lowsteps);
        this_scbins(i).meanhr_lowsteps = mean(this_scbins(i).hr_lowsteps);
        this_scbins(i).meansc_medsteps = mean(this_scbins(i).sc_medsteps);
        this_scbins(i).meanhr_medsteps = mean(this_scbins(i).hr_medsteps);
        this_scbins(i).meansc_pursteps = mean(this_scbins(i).sc_pursteps);
        this_scbins(i).meanhr_pursteps = mean(this_scbins(i).hr_pursteps);
        this_scbins(i).meansc_slowwalk = mean(this_scbins(i).sc_slowwalk);
        this_scbins(i).meanhr_slowwalk = mean(this_scbins(i).hr_slowwalk);
        this_scbins(i).meansc_medwalk = mean(this_scbins(i).sc_medwalk);
        this_scbins(i).meanhr_medwalk = mean(this_scbins(i).hr_medwalk);
        this_scbins(i).meansc_briskwalk = mean(this_scbins(i).sc_briskwalk);
        this_scbins(i).meanhr_briskwalk = mean(this_scbins(i).hr_briskwalk);
        this_scbins(i).meansc_fastamb = mean(this_scbins(i).sc_fastamb);
        this_scbins(i).meanhr_fastamb = mean(this_scbins(i).hr_fastamb);
        this_scbins(i).meansc_run = mean(this_scbins(i).sc_run);
        this_scbins(i).meanhr_run = mean(this_scbins(i).hr_run);
        
    else
        
        sc_lowsteps = []; hr_lowsteps = [];
        sc_medsteps = []; hr_medsteps = [];
        sc_pursteps = []; hr_pursteps = [];
        sc_slowwalk = []; hr_slowwalk = [];
        sc_medwalk = []; hr_medwalk = [];
        sc_briskwalk = []; hr_briskwalk = [];
        sc_fastamb = []; hr_fastamb = [];
        sc_run = []; hr_run = [];
        
        this_scbins(i).sc_lowsteps = []; this_scbins(i).hr_lowsteps = [];
        this_scbins(i).sc_medsteps = []; this_scbins(i).hr_medsteps = [];
        this_scbins(i).sc_pursteps = []; this_scbins(i).hr_pursteps = [];
        this_scbins(i).sc_slowwalk = []; this_scbins(i).hr_slowwalk = [];
        this_scbins(i).sc_medwalk = []; this_scbins(i).hr_medwalk = [];
        this_scbins(i).sc_briskwalk = []; this_scbins(i).hr_briskwalk = [];
        this_scbins(i).sc_fastamb = []; this_scbins(i).hr_fastamb = [];
        this_scbins(i).sc_run = []; this_scbins(i).hr_run = [];
            
        for j = ((this_scbins(i-1).cum_entries)+1):this_scbins(i).cum_entries
            lowsteps = find(this_weekly_raw(j).SC > 0 & this_weekly_raw(j).SC < 20);
            medsteps = find(this_weekly_raw(j).SC > 19 & this_weekly_raw(j).SC < 40);
            pursteps = find(this_weekly_raw(j).SC > 39 & this_weekly_raw(j).SC < 60);
            slowwalk = find(this_weekly_raw(j).SC > 59 & this_weekly_raw(j).SC < 80);
            medwalk = find(this_weekly_raw(j).SC > 79 & this_weekly_raw(j).SC < 100);
            briskwalk = find(this_weekly_raw(j).SC > 99 & this_weekly_raw(j).SC < 120);
            fastamb = find(this_weekly_raw(j).SC > 119 & this_weekly_raw(j).SC < 150);
            run = find(this_weekly_raw(j).SC > 149 & this_weekly_raw(j).SC < 250);
            
            if length(lowsteps) > screq - 1
                for k = 1:length(lowsteps)
                    entry = lowsteps(k);
                    sc_lowsteps = [sc_lowsteps; this_weekly_raw(j).SC(entry)];
                    hr_lowsteps = [hr_lowsteps; this_weekly_raw(j).HR(entry)];
                end
                this_scbins(i).sc_lowsteps = [this_scbins(i).sc_lowsteps; sc_lowsteps];
                this_scbins(i).hr_lowsteps = [this_scbins(i).hr_lowsteps; hr_lowsteps];
            else
                this_scbins(i).meansc_lowsteps = NaN;
                this_scbins(i).meanhr_lowsteps = NaN;
            end
            
            if length(medsteps) > screq - 1
                for k = 1:length(medsteps)
                    entry = medsteps(k);
                    sc_medsteps = [sc_medsteps; this_weekly_raw(j).SC(entry)];
                    hr_medsteps = [hr_medsteps; this_weekly_raw(j).HR(entry)];
                end
                this_scbins(i).sc_medsteps = [this_scbins(i).sc_medsteps; sc_medsteps];
                this_scbins(i).hr_medsteps = [this_scbins(i).hr_medsteps; hr_medsteps];
            else
                this_scbins(i).meansc_medsteps = NaN;
                this_scbins(i).meanhr_medsteps = NaN;
            end
            
            if length(pursteps) > screq - 1
                for k = 1:length(pursteps)
                    entry = pursteps(k);
                    sc_pursteps = [sc_pursteps; this_weekly_raw(j).SC(entry)];
                    hr_pursteps = [hr_pursteps; this_weekly_raw(j).HR(entry)];
                end
                this_scbins(i).sc_pursteps = [this_scbins(i).sc_pursteps; sc_pursteps];
                this_scbins(i).hr_pursteps = [this_scbins(i).hr_pursteps; hr_pursteps];
            else
                this_scbins(i).meansc_pursteps = NaN;
                this_scbins(i).meanhr_pursteps = NaN;       
            end
            
            if length(slowwalk) > screq - 1
                for k = 1:length(slowwalk)
                    entry = slowwalk(k);
                    sc_slowwalk = [sc_slowwalk; this_weekly_raw(j).SC(entry)];
                    hr_slowwalk = [hr_slowwalk; this_weekly_raw(j).HR(entry)];
                end
                this_scbins(i).sc_slowwalk = [this_scbins(i).sc_slowwalk; sc_slowwalk];
                this_scbins(i).hr_slowwalk = [this_scbins(i).hr_slowwalk; hr_slowwalk];
            else
                this_scbins(i).meansc_slowwalk = NaN;
                this_scbins(i).meanhr_slowwalk = NaN;
            end
            
            if length(medwalk) > screq - 1
                for k = 1:length(medwalk)
                    entry = medwalk(k);
                    sc_medwalk = [sc_medwalk; this_weekly_raw(j).SC(entry)];
                    hr_medwalk = [hr_medwalk; this_weekly_raw(j).HR(entry)];
                end
                this_scbins(i).sc_medwalk = [this_scbins(i).sc_medwalk; sc_medwalk];
                this_scbins(i).hr_medwalk = [this_scbins(i).hr_medwalk; hr_medwalk];
            else
                this_scbins(i).meansc_medwalk = NaN;
                this_scbins(i).meanhr_medwalk = NaN;
            end
            
            if length(briskwalk) > screq - 1
                for k = 1:length(briskwalk)
                    entry = briskwalk(k);
                    sc_briskwalk = [sc_briskwalk; this_weekly_raw(j).SC(entry)];
                    hr_briskwalk = [hr_briskwalk; this_weekly_raw(j).HR(entry)];
                end
                this_scbins(i).sc_briskwalk = [this_scbins(i).sc_briskwalk; sc_briskwalk];
                this_scbins(i).hr_briskwalk = [this_scbins(i).hr_briskwalk; hr_briskwalk];
            else
                this_scbins(i).meansc_briskwalk = NaN;
                this_scbins(i).meanhr_briskwalk = NaN;
            end
            
            if length(fastamb) > screq - 1
                for k = 1:length(fastamb)
                    entry = fastamb(k);
                    sc_fastamb = [sc_fastamb; this_weekly_raw(j).SC(entry)];
                    hr_fastamb = [hr_fastamb; this_weekly_raw(j).HR(entry)];
                end
                this_scbins(i).sc_fastamb = [this_scbins(i).sc_fastamb; sc_fastamb];
                this_scbins(i).hr_fastamb = [this_scbins(i).hr_fastamb; hr_fastamb];
            else
                this_scbins(i).meansc_fastamb = NaN;
                this_scbins(i).meanhr_fastamb = NaN;
            end
            
            if length(run) > screq - 1
                for k = 1:length(run)
                    entry = run(k);
                    sc_run = [sc_run; this_weekly_raw(j).SC(entry)];
                    hr_run = [hr_run; this_weekly_raw(j).HR(entry)];
                end
                this_scbins(i).sc_run = [this_scbins(i).sc_run; sc_run];
                this_scbins(i).hr_run = [this_scbins(i).hr_run; hr_run];
            else
                this_scbins(i).meansc_run = NaN;
                this_scbins(i).meanhr_run = NaN;
            end
            
        end
        
        this_scbins(i).meansc_lowsteps = mean(this_scbins(i).sc_lowsteps);
        this_scbins(i).meanhr_lowsteps = mean(this_scbins(i).hr_lowsteps);
        this_scbins(i).meansc_medsteps = mean(this_scbins(i).sc_medsteps);
        this_scbins(i).meanhr_medsteps = mean(this_scbins(i).hr_medsteps);
        this_scbins(i).meansc_pursteps = mean(this_scbins(i).sc_pursteps);
        this_scbins(i).meanhr_pursteps = mean(this_scbins(i).hr_pursteps);
        this_scbins(i).meansc_slowwalk = mean(this_scbins(i).sc_slowwalk);
        this_scbins(i).meanhr_slowwalk = mean(this_scbins(i).hr_slowwalk);
        this_scbins(i).meansc_medwalk = mean(this_scbins(i).sc_medwalk);
        this_scbins(i).meanhr_medwalk = mean(this_scbins(i).hr_medwalk);
        this_scbins(i).meansc_briskwalk = mean(this_scbins(i).sc_briskwalk);
        this_scbins(i).meanhr_briskwalk = mean(this_scbins(i).hr_briskwalk);
        this_scbins(i).meansc_fastamb = mean(this_scbins(i).sc_fastamb);
        this_scbins(i).meanhr_fastamb = mean(this_scbins(i).hr_fastamb);
        this_scbins(i).meansc_run = mean(this_scbins(i).sc_run);
        this_scbins(i).meanhr_run = mean(this_scbins(i).hr_run);
            
    end
end
    
%     this_scbins(i).sdhr_lowsteps = std(hr_lowsteps);
%     this_scbins(i).sdhr_medsteps = std(hr_medsteps);
%     this_scbins(i).sdhr_pursteps = std(hr_pursteps);
%     this_scbins(i).sdhr_slowwalk = std(hr_slowwalk);
%     this_scbins(i).sdhr_medwalk = std(hr_medwalk);
%     this_scbins(i).sdhr_briskwalk = std(hr_briskwalk);
%     this_scbins(i).sdhr_fastamb = std(hr_fastamb);
%     this_scbins(i).sdhr_run = std(hr_run);

scbins_PC = this_scbins;


%% Graph average HR in different SC bins for PC Patients

allPC_scbins = zeros(length(PC_indices),8); % Change based on # bins.
this_scbins = scbins_PC;

figure(2)
allxdata = [];
allydata = [];
hold on
for i = 1:length(PC_indices)
    xdata = [this_scbins(i).meansc_lowsteps, this_scbins(i).meansc_medsteps,...
        this_scbins(i).meansc_pursteps, this_scbins(i).meansc_slowwalk, ...
        this_scbins(i).meansc_medwalk, this_scbins(i).meansc_briskwalk, ...
        this_scbins(i).meansc_fastamb, this_scbins(i).meansc_run];
    ydata = [this_scbins(i).meanhr_lowsteps, this_scbins(i).meanhr_medsteps,...
        this_scbins(i).meanhr_pursteps, this_scbins(i).meanhr_slowwalk, ...
        this_scbins(i).meanhr_medwalk, this_scbins(i).meanhr_briskwalk, ...
        this_scbins(i).meanhr_fastamb, this_scbins(i).meanhr_run];
    allPC_scbins(i,:) = ydata;
    for j = 1:length(xdata)
        nanind = find(isnan(xdata));
        if isempty(nanind)
            xdata = xdata;
        else
            for k = 1:length(nanind)
                xdata(nanind(length(nanind)-k+1)) = [];
                ydata(nanind(length(nanind)-k+1)) = [];
            end
        end
    end
    allxdata = [allxdata, xdata];
    allydata = [allydata, ydata];
    plot(xdata,ydata,'-o','LineWidth',2)
%     err_std = [this_scbins(i).sdhr_lowsteps, this_scbins(i).sdhr_medwalk, ...
%         this_scbins(i).sdhr_pursteps, this_scbins(i).sdhr_slowwalk, ...
%         this_scbins(i).sdhr_medwalk, this_scbins(i).sdhr_briskwalk, ...
%         this_scbins(i).sdhr_fastamb, this_scbins(i).sdhr_run];
%     errorbar(xdata, ydata, err_std, '-o')
    p = polyfit(xdata,ydata, 1);
    this_scbins(i).slope = p(1);
    this_scbins(i).intercept = p(2);
end

% % Old version: using polyfit on ALL data (no differentiation b/t subjects)
% p = polyfit(allxdata,allydata, 1);
% lin_polyfit_vals(1,1) = p(1);
% lin_polyfit_vals(1,2) = p(2);
% x1 = linspace(0,max(allxdata));
% y1 = polyval(p, x1);
% plot(x1,y1,'Color','k','LineWidth',2)
% TextLocation(['y = ' num2str(p(1)), 'x + ' num2str(p(2))],'Location','southeast');

% New version: taking average slope & intercept to plot LSF line
slope = [];
int = [];
for i = 1:length(this_scbins)
    slope = [slope this_scbins(i).slope];
    int = [int this_scbins(i).intercept];
end
lin_polyfit_vals(2,1) = mean(slope);
lin_polyfit_vals(2,2) = mean(int);
x2 = linspace(0,max(allxdata));
y2 = @(x2) (mean(slope)*x2 + mean(int));
fplot(y2,'-k','LineWidth',2)
TextLocation(['y = ' num2str(mean(slope)),... 
    'x + ' num2str(mean(int))],'Location','southeast');

title('Mean HR at Various SC Thresholds for PC Patients')
xlabel('SPM')
ylabel('Mean HR (BPM)')
legend('JHPC1','JHPC9','Location','best')
hold off

scbins_PC = this_scbins;

%% Determine average HR in different SC bins for COPD Patients

% Calculate average HR for different SC bins for each COPD patient.
scbins_COPD = struct([]);
this_scbins = scbins_COPD;
this_weekly_raw = weekly_raw_COPD;

for i = 1:length(COPD_indices)
    this_scbins(i).COPD_indices = COPD_indices(i);
    this_scbins(i).num_entries = repscat_COPD(i,2);
    this_scbins(i).cum_entries = repscat_COPD(i,3);
end

for i = 1:length(COPD_indices)
    if i == 1
        
        sc_lowsteps = []; hr_lowsteps = [];
        sc_medsteps = []; hr_medsteps = [];
        sc_pursteps = []; hr_pursteps = [];
        sc_slowwalk = []; hr_slowwalk = [];
        sc_medwalk = []; hr_medwalk = [];
        sc_briskwalk = []; hr_briskwalk = [];
        sc_fastamb = []; hr_fastamb = [];
        sc_run = []; hr_run = [];
        
        this_scbins(i).sc_lowsteps = []; this_scbins(i).hr_lowsteps = [];
        this_scbins(i).sc_medsteps = []; this_scbins(i).hr_medsteps = [];
        this_scbins(i).sc_pursteps = []; this_scbins(i).hr_pursteps = [];
        this_scbins(i).sc_slowwalk = []; this_scbins(i).hr_slowwalk = [];
        this_scbins(i).sc_medwalk = []; this_scbins(i).hr_medwalk = [];
        this_scbins(i).sc_briskwalk = []; this_scbins(i).hr_briskwalk = [];
        this_scbins(i).sc_fastamb = []; this_scbins(i).hr_fastamb = [];
        this_scbins(i).sc_run = []; this_scbins(i).hr_run = [];
            
        for j = 1:this_scbins(i).cum_entries
            lowsteps = find(this_weekly_raw(j).SC > 0 & this_weekly_raw(j).SC < 20);
            medsteps = find(this_weekly_raw(j).SC > 19 & this_weekly_raw(j).SC < 40);
            pursteps = find(this_weekly_raw(j).SC > 39 & this_weekly_raw(j).SC < 60);
            slowwalk = find(this_weekly_raw(j).SC > 59 & this_weekly_raw(j).SC < 80);
            medwalk = find(this_weekly_raw(j).SC > 79 & this_weekly_raw(j).SC < 100);
            briskwalk = find(this_weekly_raw(j).SC > 99 & this_weekly_raw(j).SC < 120);
            fastamb = find(this_weekly_raw(j).SC > 119 & this_weekly_raw(j).SC < 150);
            run = find(this_weekly_raw(j).SC > 149 & this_weekly_raw(j).SC < 250);
            
            if length(lowsteps) > screq - 1
                for k = 1:length(lowsteps)
                    entry = lowsteps(k);
                    sc_lowsteps = [sc_lowsteps; this_weekly_raw(j).SC(entry)];
                    hr_lowsteps = [hr_lowsteps; this_weekly_raw(j).HR(entry)];
                end
                this_scbins(i).sc_lowsteps = [this_scbins(i).sc_lowsteps; sc_lowsteps];
                this_scbins(i).hr_lowsteps = [this_scbins(i).hr_lowsteps; hr_lowsteps];
            else
                this_scbins(i).meansc_lowsteps = NaN;
                this_scbins(i).meanhr_lowsteps = NaN;
            end
            
            if length(medsteps) > screq - 1
                for k = 1:length(medsteps)
                    entry = medsteps(k);
                    sc_medsteps = [sc_medsteps; this_weekly_raw(j).SC(entry)];
                    hr_medsteps = [hr_medsteps; this_weekly_raw(j).HR(entry)];
                end
                this_scbins(i).sc_medsteps = [this_scbins(i).sc_medsteps; sc_medsteps];
                this_scbins(i).hr_medsteps = [this_scbins(i).hr_medsteps; hr_medsteps];
            else
                this_scbins(i).meansc_medsteps = NaN;
                this_scbins(i).meanhr_medsteps = NaN;
            end
            
            if length(pursteps) > screq - 1
                for k = 1:length(pursteps)
                    entry = pursteps(k);
                    sc_pursteps = [sc_pursteps; this_weekly_raw(j).SC(entry)];
                    hr_pursteps = [hr_pursteps; this_weekly_raw(j).HR(entry)];
                end
                this_scbins(i).sc_pursteps = [this_scbins(i).sc_pursteps; sc_pursteps];
                this_scbins(i).hr_pursteps = [this_scbins(i).hr_pursteps; hr_pursteps];
            else
                this_scbins(i).meansc_pursteps = NaN;
                this_scbins(i).meanhr_pursteps = NaN;       
            end
            
            if length(slowwalk) > screq - 1
                for k = 1:length(slowwalk)
                    entry = slowwalk(k);
                    sc_slowwalk = [sc_slowwalk; this_weekly_raw(j).SC(entry)];
                    hr_slowwalk = [hr_slowwalk; this_weekly_raw(j).HR(entry)];
                end
                this_scbins(i).sc_slowwalk = [this_scbins(i).sc_slowwalk; sc_slowwalk];
                this_scbins(i).hr_slowwalk = [this_scbins(i).hr_slowwalk; hr_slowwalk];
            else
                this_scbins(i).meansc_slowwalk = NaN;
                this_scbins(i).meanhr_slowwalk = NaN;
            end
            
            if length(medwalk) > screq - 1
                for k = 1:length(medwalk)
                    entry = medwalk(k);
                    sc_medwalk = [sc_medwalk; this_weekly_raw(j).SC(entry)];
                    hr_medwalk = [hr_medwalk; this_weekly_raw(j).HR(entry)];
                end
                this_scbins(i).sc_medwalk = [this_scbins(i).sc_medwalk; sc_medwalk];
                this_scbins(i).hr_medwalk = [this_scbins(i).hr_medwalk; hr_medwalk];
            else
                this_scbins(i).meansc_medwalk = NaN;
                this_scbins(i).meanhr_medwalk = NaN;
            end
            
            if length(briskwalk) > screq - 1
                for k = 1:length(briskwalk)
                    entry = briskwalk(k);
                    sc_briskwalk = [sc_briskwalk; this_weekly_raw(j).SC(entry)];
                    hr_briskwalk = [hr_briskwalk; this_weekly_raw(j).HR(entry)];
                end
                this_scbins(i).sc_briskwalk = [this_scbins(i).sc_briskwalk; sc_briskwalk];
                this_scbins(i).hr_briskwalk = [this_scbins(i).hr_briskwalk; hr_briskwalk];
            else
                this_scbins(i).meansc_briskwalk = NaN;
                this_scbins(i).meanhr_briskwalk = NaN;
            end
            
            if length(fastamb) > screq - 1
                for k = 1:length(fastamb)
                    entry = fastamb(k);
                    sc_fastamb = [sc_fastamb; this_weekly_raw(j).SC(entry)];
                    hr_fastamb = [hr_fastamb; this_weekly_raw(j).HR(entry)];
                end
                this_scbins(i).sc_fastamb = [this_scbins(i).sc_fastamb; sc_fastamb];
                this_scbins(i).hr_fastamb = [this_scbins(i).hr_fastamb; hr_fastamb];
            else
                this_scbins(i).meansc_fastamb = NaN;
                this_scbins(i).meanhr_fastamb = NaN;
            end
            
            if length(run) > screq - 1
                for k = 1:length(run)
                    entry = run(k);
                    sc_run = [sc_run; this_weekly_raw(j).SC(entry)];
                    hr_run = [hr_run; this_weekly_raw(j).HR(entry)];
                end
                this_scbins(i).sc_run = [this_scbins(i).sc_run; sc_run];
                this_scbins(i).hr_run = [this_scbins(i).hr_run; hr_run];
            else
                this_scbins(i).meansc_run = NaN;
                this_scbins(i).meanhr_run = NaN;
            end
            
        end
        
        this_scbins(i).meansc_lowsteps = mean(this_scbins(i).sc_lowsteps);
        this_scbins(i).meanhr_lowsteps = mean(this_scbins(i).hr_lowsteps);
        this_scbins(i).meansc_medsteps = mean(this_scbins(i).sc_medsteps);
        this_scbins(i).meanhr_medsteps = mean(this_scbins(i).hr_medsteps);
        this_scbins(i).meansc_pursteps = mean(this_scbins(i).sc_pursteps);
        this_scbins(i).meanhr_pursteps = mean(this_scbins(i).hr_pursteps);
        this_scbins(i).meansc_slowwalk = mean(this_scbins(i).sc_slowwalk);
        this_scbins(i).meanhr_slowwalk = mean(this_scbins(i).hr_slowwalk);
        this_scbins(i).meansc_medwalk = mean(this_scbins(i).sc_medwalk);
        this_scbins(i).meanhr_medwalk = mean(this_scbins(i).hr_medwalk);
        this_scbins(i).meansc_briskwalk = mean(this_scbins(i).sc_briskwalk);
        this_scbins(i).meanhr_briskwalk = mean(this_scbins(i).hr_briskwalk);
        this_scbins(i).meansc_fastamb = mean(this_scbins(i).sc_fastamb);
        this_scbins(i).meanhr_fastamb = mean(this_scbins(i).hr_fastamb);
        this_scbins(i).meansc_run = mean(this_scbins(i).sc_run);
        this_scbins(i).meanhr_run = mean(this_scbins(i).hr_run);
        
    else
        
        sc_lowsteps = []; hr_lowsteps = [];
        sc_medsteps = []; hr_medsteps = [];
        sc_pursteps = []; hr_pursteps = [];
        sc_slowwalk = []; hr_slowwalk = [];
        sc_medwalk = []; hr_medwalk = [];
        sc_briskwalk = []; hr_briskwalk = [];
        sc_fastamb = []; hr_fastamb = [];
        sc_run = []; hr_run = [];
        
        this_scbins(i).sc_lowsteps = []; this_scbins(i).hr_lowsteps = [];
        this_scbins(i).sc_medsteps = []; this_scbins(i).hr_medsteps = [];
        this_scbins(i).sc_pursteps = []; this_scbins(i).hr_pursteps = [];
        this_scbins(i).sc_slowwalk = []; this_scbins(i).hr_slowwalk = [];
        this_scbins(i).sc_medwalk = []; this_scbins(i).hr_medwalk = [];
        this_scbins(i).sc_briskwalk = []; this_scbins(i).hr_briskwalk = [];
        this_scbins(i).sc_fastamb = []; this_scbins(i).hr_fastamb = [];
        this_scbins(i).sc_run = []; this_scbins(i).hr_run = [];
            
        for j = ((this_scbins(i-1).cum_entries)+1):this_scbins(i).cum_entries
            lowsteps = find(this_weekly_raw(j).SC > 0 & this_weekly_raw(j).SC < 20);
            medsteps = find(this_weekly_raw(j).SC > 19 & this_weekly_raw(j).SC < 40);
            pursteps = find(this_weekly_raw(j).SC > 39 & this_weekly_raw(j).SC < 60);
            slowwalk = find(this_weekly_raw(j).SC > 59 & this_weekly_raw(j).SC < 80);
            medwalk = find(this_weekly_raw(j).SC > 79 & this_weekly_raw(j).SC < 100);
            briskwalk = find(this_weekly_raw(j).SC > 99 & this_weekly_raw(j).SC < 120);
            fastamb = find(this_weekly_raw(j).SC > 119 & this_weekly_raw(j).SC < 150);
            run = find(this_weekly_raw(j).SC > 149 & this_weekly_raw(j).SC < 250);
            
            if length(lowsteps) > screq - 1
                for k = 1:length(lowsteps)
                    entry = lowsteps(k);
                    sc_lowsteps = [sc_lowsteps; this_weekly_raw(j).SC(entry)];
                    hr_lowsteps = [hr_lowsteps; this_weekly_raw(j).HR(entry)];
                end
                this_scbins(i).sc_lowsteps = [this_scbins(i).sc_lowsteps; sc_lowsteps];
                this_scbins(i).hr_lowsteps = [this_scbins(i).hr_lowsteps; hr_lowsteps];
            else
                this_scbins(i).meansc_lowsteps = NaN;
                this_scbins(i).meanhr_lowsteps = NaN;
            end
            
            if length(medsteps) > screq - 1
                for k = 1:length(medsteps)
                    entry = medsteps(k);
                    sc_medsteps = [sc_medsteps; this_weekly_raw(j).SC(entry)];
                    hr_medsteps = [hr_medsteps; this_weekly_raw(j).HR(entry)];
                end
                this_scbins(i).sc_medsteps = [this_scbins(i).sc_medsteps; sc_medsteps];
                this_scbins(i).hr_medsteps = [this_scbins(i).hr_medsteps; hr_medsteps];
            else
                this_scbins(i).meansc_medsteps = NaN;
                this_scbins(i).meanhr_medsteps = NaN;
            end
            
            if length(pursteps) > screq - 1
                for k = 1:length(pursteps)
                    entry = pursteps(k);
                    sc_pursteps = [sc_pursteps; this_weekly_raw(j).SC(entry)];
                    hr_pursteps = [hr_pursteps; this_weekly_raw(j).HR(entry)];
                end
                this_scbins(i).sc_pursteps = [this_scbins(i).sc_pursteps; sc_pursteps];
                this_scbins(i).hr_pursteps = [this_scbins(i).hr_pursteps; hr_pursteps];
            else
                this_scbins(i).meansc_pursteps = NaN;
                this_scbins(i).meanhr_pursteps = NaN;       
            end
            
            if length(slowwalk) > screq - 1
                for k = 1:length(slowwalk)
                    entry = slowwalk(k);
                    sc_slowwalk = [sc_slowwalk; this_weekly_raw(j).SC(entry)];
                    hr_slowwalk = [hr_slowwalk; this_weekly_raw(j).HR(entry)];
                end
                this_scbins(i).sc_slowwalk = [this_scbins(i).sc_slowwalk; sc_slowwalk];
                this_scbins(i).hr_slowwalk = [this_scbins(i).hr_slowwalk; hr_slowwalk];
            else
                this_scbins(i).meansc_slowwalk = NaN;
                this_scbins(i).meanhr_slowwalk = NaN;
            end
            
            if length(medwalk) > screq - 1
                for k = 1:length(medwalk)
                    entry = medwalk(k);
                    sc_medwalk = [sc_medwalk; this_weekly_raw(j).SC(entry)];
                    hr_medwalk = [hr_medwalk; this_weekly_raw(j).HR(entry)];
                end
                this_scbins(i).sc_medwalk = [this_scbins(i).sc_medwalk; sc_medwalk];
                this_scbins(i).hr_medwalk = [this_scbins(i).hr_medwalk; hr_medwalk];
            else
                this_scbins(i).meansc_medwalk = NaN;
                this_scbins(i).meanhr_medwalk = NaN;
            end
            
            if length(briskwalk) > screq - 1
                for k = 1:length(briskwalk)
                    entry = briskwalk(k);
                    sc_briskwalk = [sc_briskwalk; this_weekly_raw(j).SC(entry)];
                    hr_briskwalk = [hr_briskwalk; this_weekly_raw(j).HR(entry)];
                end
                this_scbins(i).sc_briskwalk = [this_scbins(i).sc_briskwalk; sc_briskwalk];
                this_scbins(i).hr_briskwalk = [this_scbins(i).hr_briskwalk; hr_briskwalk];
            else
                this_scbins(i).meansc_briskwalk = NaN;
                this_scbins(i).meanhr_briskwalk = NaN;
            end
            
            if length(fastamb) > screq - 1
                for k = 1:length(fastamb)
                    entry = fastamb(k);
                    sc_fastamb = [sc_fastamb; this_weekly_raw(j).SC(entry)];
                    hr_fastamb = [hr_fastamb; this_weekly_raw(j).HR(entry)];
                end
                this_scbins(i).sc_fastamb = [this_scbins(i).sc_fastamb; sc_fastamb];
                this_scbins(i).hr_fastamb = [this_scbins(i).hr_fastamb; hr_fastamb];
            else
                this_scbins(i).meansc_fastamb = NaN;
                this_scbins(i).meanhr_fastamb = NaN;
            end
            
            if length(run) > screq - 1
                for k = 1:length(run)
                    entry = run(k);
                    sc_run = [sc_run; this_weekly_raw(j).SC(entry)];
                    hr_run = [hr_run; this_weekly_raw(j).HR(entry)];
                end
                this_scbins(i).sc_run = [this_scbins(i).sc_run; sc_run];
                this_scbins(i).hr_run = [this_scbins(i).hr_run; hr_run];
            else
                this_scbins(i).meansc_run = NaN;
                this_scbins(i).meanhr_run = NaN;
            end
            
        end
        
        this_scbins(i).meansc_lowsteps = mean(this_scbins(i).sc_lowsteps);
        this_scbins(i).meanhr_lowsteps = mean(this_scbins(i).hr_lowsteps);
        this_scbins(i).meansc_medsteps = mean(this_scbins(i).sc_medsteps);
        this_scbins(i).meanhr_medsteps = mean(this_scbins(i).hr_medsteps);
        this_scbins(i).meansc_pursteps = mean(this_scbins(i).sc_pursteps);
        this_scbins(i).meanhr_pursteps = mean(this_scbins(i).hr_pursteps);
        this_scbins(i).meansc_slowwalk = mean(this_scbins(i).sc_slowwalk);
        this_scbins(i).meanhr_slowwalk = mean(this_scbins(i).hr_slowwalk);
        this_scbins(i).meansc_medwalk = mean(this_scbins(i).sc_medwalk);
        this_scbins(i).meanhr_medwalk = mean(this_scbins(i).hr_medwalk);
        this_scbins(i).meansc_briskwalk = mean(this_scbins(i).sc_briskwalk);
        this_scbins(i).meanhr_briskwalk = mean(this_scbins(i).hr_briskwalk);
        this_scbins(i).meansc_fastamb = mean(this_scbins(i).sc_fastamb);
        this_scbins(i).meanhr_fastamb = mean(this_scbins(i).hr_fastamb);
        this_scbins(i).meansc_run = mean(this_scbins(i).sc_run);
        this_scbins(i).meanhr_run = mean(this_scbins(i).hr_run);
            
    end
end
    
%     this_scbins(i).sdhr_lowsteps = std(hr_lowsteps);
%     this_scbins(i).sdhr_medsteps = std(hr_medsteps);
%     this_scbins(i).sdhr_pursteps = std(hr_pursteps);
%     this_scbins(i).sdhr_slowwalk = std(hr_slowwalk);
%     this_scbins(i).sdhr_medwalk = std(hr_medwalk);
%     this_scbins(i).sdhr_briskwalk = std(hr_briskwalk);
%     this_scbins(i).sdhr_fastamb = std(hr_fastamb);
%     this_scbins(i).sdhr_run = std(hr_run);

scbins_COPD = this_scbins;


%% Graph average HR in different SC bins for COPD Patients

allCOPD_scbins = zeros(length(COPD_indices),8); % Change based on # bins.
this_scbins = scbins_COPD;
 
figure(3)
allxdata = [];
allydata = [];
hold on
for i = 1:length(COPD_indices)
    xdata = [this_scbins(i).meansc_lowsteps, this_scbins(i).meansc_medsteps,...
        this_scbins(i).meansc_pursteps, this_scbins(i).meansc_slowwalk, ...
        this_scbins(i).meansc_medwalk, this_scbins(i).meansc_briskwalk, ...
        this_scbins(i).meansc_fastamb, this_scbins(i).meansc_run];
    ydata = [this_scbins(i).meanhr_lowsteps, this_scbins(i).meanhr_medsteps,...
        this_scbins(i).meanhr_pursteps, this_scbins(i).meanhr_slowwalk, ...
        this_scbins(i).meanhr_medwalk, this_scbins(i).meanhr_briskwalk, ...
        this_scbins(i).meanhr_fastamb, this_scbins(i).meanhr_run];
    allCOPD_scbins(i,:) = ydata;
    for j = 1:length(xdata)
        nanind = find(isnan(xdata));
        if isempty(nanind)
            xdata = xdata;
        else
            for k = 1:length(nanind)
                xdata(nanind(length(nanind)-k+1)) = [];
                ydata(nanind(length(nanind)-k+1)) = [];
            end
        end
    end
    allxdata = [allxdata, xdata];
    allydata = [allydata, ydata];
    plot(xdata,ydata,'-o','LineWidth',2)
%     err_std = [this_scbins(i).sdhr_lowsteps, this_scbins(i).sdhr_medwalk, ...
%         this_scbins(i).sdhr_pursteps, this_scbins(i).sdhr_slowwalk, ...
%         this_scbins(i).sdhr_medwalk, this_scbins(i).sdhr_briskwalk, ...
%         this_scbins(i).sdhr_fastamb, this_scbins(i).sdhr_run];
%     errorbar(xdata, ydata, err_std, '-o')
    p = polyfit(xdata,ydata, 1);
    this_scbins(i).slope = p(1);
    this_scbins(i).intercept = p(2);
end

% % Old version: using polyfit on ALL data (no differentiation b/t subjects)
% p = polyfit(allxdata,allydata, 1);
% lin_polyfit_vals(1,1) = p(1);
% lin_polyfit_vals(1,2) = p(2);
% x1 = linspace(0,max(allxdata));
% y1 = polyval(p, x1);
% plot(x1,y1,'Color','k','LineWidth',2)
% TextLocation(['y = ' num2str(p(1)), 'x + ' num2str(p(2))],'Location','southeast');

% New version: taking average slope & intercept to plot LSF line
slope = [];
int = [];
for i = 1:length(this_scbins)
    slope = [slope this_scbins(i).slope];
    int = [int this_scbins(i).intercept];
end
lin_polyfit_vals(3,1) = mean(slope);
lin_polyfit_vals(3,2) = mean(int);
x2 = linspace(0,max(allxdata));
y2 = @(x2) (mean(slope)*x2 + mean(int));
fplot(y2,'-k','LineWidth',2)
TextLocation(['y = ' num2str(mean(slope)),... 
    'x + ' num2str(mean(int))],'Location','southeast');

title('Mean HR at Various SC Thresholds for COPD Patients')
xlabel('SPM')
ylabel('Mean HR (BPM)')
legend('COPD2','COPD3','COPD4','COPD5','COPD6','Location','best')
hold off

scbins_COPD = this_scbins;

%% Determine average HR in different SC bins for PAH Patients

% Calculate average HR for different SC bins for each PAH patient.
scbins_PAH = struct([]);
this_scbins = scbins_PAH;
this_weekly_raw = weekly_raw_PAH;

for i = 1:length(PAH_indices)
    this_scbins(i).PAH_indices = PAH_indices(i);
    this_scbins(i).num_entries = repscat_PAH(i,2);
    this_scbins(i).cum_entries = repscat_PAH(i,3);
end
 
for i = 1:length(PAH_indices)
    if i == 1
        
        sc_lowsteps = []; hr_lowsteps = [];
        sc_medsteps = []; hr_medsteps = [];
        sc_pursteps = []; hr_pursteps = [];
        sc_slowwalk = []; hr_slowwalk = [];
        sc_medwalk = []; hr_medwalk = [];
        sc_briskwalk = []; hr_briskwalk = [];
        sc_fastamb = []; hr_fastamb = [];
        sc_run = []; hr_run = [];
        
        this_scbins(i).sc_lowsteps = []; this_scbins(i).hr_lowsteps = [];
        this_scbins(i).sc_medsteps = []; this_scbins(i).hr_medsteps = [];
        this_scbins(i).sc_pursteps = []; this_scbins(i).hr_pursteps = [];
        this_scbins(i).sc_slowwalk = []; this_scbins(i).hr_slowwalk = [];
        this_scbins(i).sc_medwalk = []; this_scbins(i).hr_medwalk = [];
        this_scbins(i).sc_briskwalk = []; this_scbins(i).hr_briskwalk = [];
        this_scbins(i).sc_fastamb = []; this_scbins(i).hr_fastamb = [];
        this_scbins(i).sc_run = []; this_scbins(i).hr_run = [];
            
        for j = 1:this_scbins(i).cum_entries
            lowsteps = find(this_weekly_raw(j).SC > 0 & this_weekly_raw(j).SC < 20);
            medsteps = find(this_weekly_raw(j).SC > 19 & this_weekly_raw(j).SC < 40);
            pursteps = find(this_weekly_raw(j).SC > 39 & this_weekly_raw(j).SC < 60);
            slowwalk = find(this_weekly_raw(j).SC > 59 & this_weekly_raw(j).SC < 80);
            medwalk = find(this_weekly_raw(j).SC > 79 & this_weekly_raw(j).SC < 100);
            briskwalk = find(this_weekly_raw(j).SC > 99 & this_weekly_raw(j).SC < 120);
            fastamb = find(this_weekly_raw(j).SC > 119 & this_weekly_raw(j).SC < 150);
            run = find(this_weekly_raw(j).SC > 149 & this_weekly_raw(j).SC < 250);
            
            if length(lowsteps) > screq - 1
                for k = 1:length(lowsteps)
                    entry = lowsteps(k);
                    sc_lowsteps = [sc_lowsteps; this_weekly_raw(j).SC(entry)];
                    hr_lowsteps = [hr_lowsteps; this_weekly_raw(j).HR(entry)];
                end
                this_scbins(i).sc_lowsteps = [this_scbins(i).sc_lowsteps; sc_lowsteps];
                this_scbins(i).hr_lowsteps = [this_scbins(i).hr_lowsteps; hr_lowsteps];
            else
                this_scbins(i).meansc_lowsteps = NaN;
                this_scbins(i).meanhr_lowsteps = NaN;
            end
            
            if length(medsteps) > screq - 1
                for k = 1:length(medsteps)
                    entry = medsteps(k);
                    sc_medsteps = [sc_medsteps; this_weekly_raw(j).SC(entry)];
                    hr_medsteps = [hr_medsteps; this_weekly_raw(j).HR(entry)];
                end
                this_scbins(i).sc_medsteps = [this_scbins(i).sc_medsteps; sc_medsteps];
                this_scbins(i).hr_medsteps = [this_scbins(i).hr_medsteps; hr_medsteps];
            else
                this_scbins(i).meansc_medsteps = NaN;
                this_scbins(i).meanhr_medsteps = NaN;
            end
            
            if length(pursteps) > screq - 1
                for k = 1:length(pursteps)
                    entry = pursteps(k);
                    sc_pursteps = [sc_pursteps; this_weekly_raw(j).SC(entry)];
                    hr_pursteps = [hr_pursteps; this_weekly_raw(j).HR(entry)];
                end
                this_scbins(i).sc_pursteps = [this_scbins(i).sc_pursteps; sc_pursteps];
                this_scbins(i).hr_pursteps = [this_scbins(i).hr_pursteps; hr_pursteps];
            else
                this_scbins(i).meansc_pursteps = NaN;
                this_scbins(i).meanhr_pursteps = NaN;       
            end
            
            if length(slowwalk) > screq - 1
                for k = 1:length(slowwalk)
                    entry = slowwalk(k);
                    sc_slowwalk = [sc_slowwalk; this_weekly_raw(j).SC(entry)];
                    hr_slowwalk = [hr_slowwalk; this_weekly_raw(j).HR(entry)];
                end
                this_scbins(i).sc_slowwalk = [this_scbins(i).sc_slowwalk; sc_slowwalk];
                this_scbins(i).hr_slowwalk = [this_scbins(i).hr_slowwalk; hr_slowwalk];
            else
                this_scbins(i).meansc_slowwalk = NaN;
                this_scbins(i).meanhr_slowwalk = NaN;
            end
            
            if length(medwalk) > screq - 1
                for k = 1:length(medwalk)
                    entry = medwalk(k);
                    sc_medwalk = [sc_medwalk; this_weekly_raw(j).SC(entry)];
                    hr_medwalk = [hr_medwalk; this_weekly_raw(j).HR(entry)];
                end
                this_scbins(i).sc_medwalk = [this_scbins(i).sc_medwalk; sc_medwalk];
                this_scbins(i).hr_medwalk = [this_scbins(i).hr_medwalk; hr_medwalk];
            else
                this_scbins(i).meansc_medwalk = NaN;
                this_scbins(i).meanhr_medwalk = NaN;
            end
            
            if length(briskwalk) > screq - 1
                for k = 1:length(briskwalk)
                    entry = briskwalk(k);
                    sc_briskwalk = [sc_briskwalk; this_weekly_raw(j).SC(entry)];
                    hr_briskwalk = [hr_briskwalk; this_weekly_raw(j).HR(entry)];
                end
                this_scbins(i).sc_briskwalk = [this_scbins(i).sc_briskwalk; sc_briskwalk];
                this_scbins(i).hr_briskwalk = [this_scbins(i).hr_briskwalk; hr_briskwalk];
            else
                this_scbins(i).meansc_briskwalk = NaN;
                this_scbins(i).meanhr_briskwalk = NaN;
            end
            
            if length(fastamb) > screq - 1
                for k = 1:length(fastamb)
                    entry = fastamb(k);
                    sc_fastamb = [sc_fastamb; this_weekly_raw(j).SC(entry)];
                    hr_fastamb = [hr_fastamb; this_weekly_raw(j).HR(entry)];
                end
                this_scbins(i).sc_fastamb = [this_scbins(i).sc_fastamb; sc_fastamb];
                this_scbins(i).hr_fastamb = [this_scbins(i).hr_fastamb; hr_fastamb];
            else
                this_scbins(i).meansc_fastamb = NaN;
                this_scbins(i).meanhr_fastamb = NaN;
            end
            
            if length(run) > screq - 1
                for k = 1:length(run)
                    entry = run(k);
                    sc_run = [sc_run; this_weekly_raw(j).SC(entry)];
                    hr_run = [hr_run; this_weekly_raw(j).HR(entry)];
                end
                this_scbins(i).sc_run = [this_scbins(i).sc_run; sc_run];
                this_scbins(i).hr_run = [this_scbins(i).hr_run; hr_run];
            else
                this_scbins(i).meansc_run = NaN;
                this_scbins(i).meanhr_run = NaN;
            end
            
        end
        
        this_scbins(i).meansc_lowsteps = mean(this_scbins(i).sc_lowsteps);
        this_scbins(i).meanhr_lowsteps = mean(this_scbins(i).hr_lowsteps);
        this_scbins(i).meansc_medsteps = mean(this_scbins(i).sc_medsteps);
        this_scbins(i).meanhr_medsteps = mean(this_scbins(i).hr_medsteps);
        this_scbins(i).meansc_pursteps = mean(this_scbins(i).sc_pursteps);
        this_scbins(i).meanhr_pursteps = mean(this_scbins(i).hr_pursteps);
        this_scbins(i).meansc_slowwalk = mean(this_scbins(i).sc_slowwalk);
        this_scbins(i).meanhr_slowwalk = mean(this_scbins(i).hr_slowwalk);
        this_scbins(i).meansc_medwalk = mean(this_scbins(i).sc_medwalk);
        this_scbins(i).meanhr_medwalk = mean(this_scbins(i).hr_medwalk);
        this_scbins(i).meansc_briskwalk = mean(this_scbins(i).sc_briskwalk);
        this_scbins(i).meanhr_briskwalk = mean(this_scbins(i).hr_briskwalk);
        this_scbins(i).meansc_fastamb = mean(this_scbins(i).sc_fastamb);
        this_scbins(i).meanhr_fastamb = mean(this_scbins(i).hr_fastamb);
        this_scbins(i).meansc_run = mean(this_scbins(i).sc_run);
        this_scbins(i).meanhr_run = mean(this_scbins(i).hr_run);
        
    else
        
        sc_lowsteps = []; hr_lowsteps = [];
        sc_medsteps = []; hr_medsteps = [];
        sc_pursteps = []; hr_pursteps = [];
        sc_slowwalk = []; hr_slowwalk = [];
        sc_medwalk = []; hr_medwalk = [];
        sc_briskwalk = []; hr_briskwalk = [];
        sc_fastamb = []; hr_fastamb = [];
        sc_run = []; hr_run = [];
        
        this_scbins(i).sc_lowsteps = []; this_scbins(i).hr_lowsteps = [];
        this_scbins(i).sc_medsteps = []; this_scbins(i).hr_medsteps = [];
        this_scbins(i).sc_pursteps = []; this_scbins(i).hr_pursteps = [];
        this_scbins(i).sc_slowwalk = []; this_scbins(i).hr_slowwalk = [];
        this_scbins(i).sc_medwalk = []; this_scbins(i).hr_medwalk = [];
        this_scbins(i).sc_briskwalk = []; this_scbins(i).hr_briskwalk = [];
        this_scbins(i).sc_fastamb = []; this_scbins(i).hr_fastamb = [];
        this_scbins(i).sc_run = []; this_scbins(i).hr_run = [];
            
        for j = ((this_scbins(i-1).cum_entries)+1):this_scbins(i).cum_entries
            lowsteps = find(this_weekly_raw(j).SC > 0 & this_weekly_raw(j).SC < 20);
            medsteps = find(this_weekly_raw(j).SC > 19 & this_weekly_raw(j).SC < 40);
            pursteps = find(this_weekly_raw(j).SC > 39 & this_weekly_raw(j).SC < 60);
            slowwalk = find(this_weekly_raw(j).SC > 59 & this_weekly_raw(j).SC < 80);
            medwalk = find(this_weekly_raw(j).SC > 79 & this_weekly_raw(j).SC < 100);
            briskwalk = find(this_weekly_raw(j).SC > 99 & this_weekly_raw(j).SC < 120);
            fastamb = find(this_weekly_raw(j).SC > 119 & this_weekly_raw(j).SC < 150);
            run = find(this_weekly_raw(j).SC > 149 & this_weekly_raw(j).SC < 250);
            
            if length(lowsteps) > screq - 1
                for k = 1:length(lowsteps)
                    entry = lowsteps(k);
                    sc_lowsteps = [sc_lowsteps; this_weekly_raw(j).SC(entry)];
                    hr_lowsteps = [hr_lowsteps; this_weekly_raw(j).HR(entry)];
                end
                this_scbins(i).sc_lowsteps = [this_scbins(i).sc_lowsteps; sc_lowsteps];
                this_scbins(i).hr_lowsteps = [this_scbins(i).hr_lowsteps; hr_lowsteps];
            else
                this_scbins(i).meansc_lowsteps = NaN;
                this_scbins(i).meanhr_lowsteps = NaN;
            end
            
            if length(medsteps) > screq - 1
                for k = 1:length(medsteps)
                    entry = medsteps(k);
                    sc_medsteps = [sc_medsteps; this_weekly_raw(j).SC(entry)];
                    hr_medsteps = [hr_medsteps; this_weekly_raw(j).HR(entry)];
                end
                this_scbins(i).sc_medsteps = [this_scbins(i).sc_medsteps; sc_medsteps];
                this_scbins(i).hr_medsteps = [this_scbins(i).hr_medsteps; hr_medsteps];
            else
                this_scbins(i).meansc_medsteps = NaN;
                this_scbins(i).meanhr_medsteps = NaN;
            end
            
            if length(pursteps) > screq - 1
                for k = 1:length(pursteps)
                    entry = pursteps(k);
                    sc_pursteps = [sc_pursteps; this_weekly_raw(j).SC(entry)];
                    hr_pursteps = [hr_pursteps; this_weekly_raw(j).HR(entry)];
                end
                this_scbins(i).sc_pursteps = [this_scbins(i).sc_pursteps; sc_pursteps];
                this_scbins(i).hr_pursteps = [this_scbins(i).hr_pursteps; hr_pursteps];
            else
                this_scbins(i).meansc_pursteps = NaN;
                this_scbins(i).meanhr_pursteps = NaN;       
            end
            
            if length(slowwalk) > screq - 1
                for k = 1:length(slowwalk)
                    entry = slowwalk(k);
                    sc_slowwalk = [sc_slowwalk; this_weekly_raw(j).SC(entry)];
                    hr_slowwalk = [hr_slowwalk; this_weekly_raw(j).HR(entry)];
                end
                this_scbins(i).sc_slowwalk = [this_scbins(i).sc_slowwalk; sc_slowwalk];
                this_scbins(i).hr_slowwalk = [this_scbins(i).hr_slowwalk; hr_slowwalk];
            else
                this_scbins(i).meansc_slowwalk = NaN;
                this_scbins(i).meanhr_slowwalk = NaN;
            end
            
            if length(medwalk) > screq - 1
                for k = 1:length(medwalk)
                    entry = medwalk(k);
                    sc_medwalk = [sc_medwalk; this_weekly_raw(j).SC(entry)];
                    hr_medwalk = [hr_medwalk; this_weekly_raw(j).HR(entry)];
                end
                this_scbins(i).sc_medwalk = [this_scbins(i).sc_medwalk; sc_medwalk];
                this_scbins(i).hr_medwalk = [this_scbins(i).hr_medwalk; hr_medwalk];
            else
                this_scbins(i).meansc_medwalk = NaN;
                this_scbins(i).meanhr_medwalk = NaN;
            end
            
            if length(briskwalk) > screq - 1
                for k = 1:length(briskwalk)
                    entry = briskwalk(k);
                    sc_briskwalk = [sc_briskwalk; this_weekly_raw(j).SC(entry)];
                    hr_briskwalk = [hr_briskwalk; this_weekly_raw(j).HR(entry)];
                end
                this_scbins(i).sc_briskwalk = [this_scbins(i).sc_briskwalk; sc_briskwalk];
                this_scbins(i).hr_briskwalk = [this_scbins(i).hr_briskwalk; hr_briskwalk];
            else
                this_scbins(i).meansc_briskwalk = NaN;
                this_scbins(i).meanhr_briskwalk = NaN;
            end
            
            if length(fastamb) > screq - 1
                for k = 1:length(fastamb)
                    entry = fastamb(k);
                    sc_fastamb = [sc_fastamb; this_weekly_raw(j).SC(entry)];
                    hr_fastamb = [hr_fastamb; this_weekly_raw(j).HR(entry)];
                end
                this_scbins(i).sc_fastamb = [this_scbins(i).sc_fastamb; sc_fastamb];
                this_scbins(i).hr_fastamb = [this_scbins(i).hr_fastamb; hr_fastamb];
            else
                this_scbins(i).meansc_fastamb = NaN;
                this_scbins(i).meanhr_fastamb = NaN;
            end
            
            if length(run) > screq - 1
                for k = 1:length(run)
                    entry = run(k);
                    sc_run = [sc_run; this_weekly_raw(j).SC(entry)];
                    hr_run = [hr_run; this_weekly_raw(j).HR(entry)];
                end
                this_scbins(i).sc_run = [this_scbins(i).sc_run; sc_run];
                this_scbins(i).hr_run = [this_scbins(i).hr_run; hr_run];
            else
                this_scbins(i).meansc_run = NaN;
                this_scbins(i).meanhr_run = NaN;
            end
            
        end
        
        this_scbins(i).meansc_lowsteps = mean(this_scbins(i).sc_lowsteps);
        this_scbins(i).meanhr_lowsteps = mean(this_scbins(i).hr_lowsteps);
        this_scbins(i).meansc_medsteps = mean(this_scbins(i).sc_medsteps);
        this_scbins(i).meanhr_medsteps = mean(this_scbins(i).hr_medsteps);
        this_scbins(i).meansc_pursteps = mean(this_scbins(i).sc_pursteps);
        this_scbins(i).meanhr_pursteps = mean(this_scbins(i).hr_pursteps);
        this_scbins(i).meansc_slowwalk = mean(this_scbins(i).sc_slowwalk);
        this_scbins(i).meanhr_slowwalk = mean(this_scbins(i).hr_slowwalk);
        this_scbins(i).meansc_medwalk = mean(this_scbins(i).sc_medwalk);
        this_scbins(i).meanhr_medwalk = mean(this_scbins(i).hr_medwalk);
        this_scbins(i).meansc_briskwalk = mean(this_scbins(i).sc_briskwalk);
        this_scbins(i).meanhr_briskwalk = mean(this_scbins(i).hr_briskwalk);
        this_scbins(i).meansc_fastamb = mean(this_scbins(i).sc_fastamb);
        this_scbins(i).meanhr_fastamb = mean(this_scbins(i).hr_fastamb);
        this_scbins(i).meansc_run = mean(this_scbins(i).sc_run);
        this_scbins(i).meanhr_run = mean(this_scbins(i).hr_run);
            
    end
end
    
%     this_scbins(i).sdhr_lowsteps = std(hr_lowsteps);
%     this_scbins(i).sdhr_medsteps = std(hr_medsteps);
%     this_scbins(i).sdhr_pursteps = std(hr_pursteps);
%     this_scbins(i).sdhr_slowwalk = std(hr_slowwalk);
%     this_scbins(i).sdhr_medwalk = std(hr_medwalk);
%     this_scbins(i).sdhr_briskwalk = std(hr_briskwalk);
%     this_scbins(i).sdhr_fastamb = std(hr_fastamb);
%     this_scbins(i).sdhr_run = std(hr_run);

scbins_PAH = this_scbins;

%% Graph average HR in different SC bins for PAH Patients

allPAH_scbins = zeros(length(PAH_indices),8); % Change based on # bins.
this_scbins = scbins_PAH;
 
figure(4)
allxdata = [];
allydata = [];
hold on
for i = 1:length(PAH_indices)
    xdata = [this_scbins(i).meansc_lowsteps, this_scbins(i).meansc_medsteps,...
        this_scbins(i).meansc_pursteps, this_scbins(i).meansc_slowwalk, ...
        this_scbins(i).meansc_medwalk, this_scbins(i).meansc_briskwalk, ...
        this_scbins(i).meansc_fastamb, this_scbins(i).meansc_run];
    ydata = [this_scbins(i).meanhr_lowsteps, this_scbins(i).meanhr_medsteps,...
        this_scbins(i).meanhr_pursteps, this_scbins(i).meanhr_slowwalk, ...
        this_scbins(i).meanhr_medwalk, this_scbins(i).meanhr_briskwalk, ...
        this_scbins(i).meanhr_fastamb, this_scbins(i).meanhr_run];
    allPAH_scbins(i,:) = ydata;
    for j = 1:length(xdata)
        nanind = find(isnan(xdata));
        if isempty(nanind)
            xdata = xdata;
        else
            for k = 1:length(nanind)
                xdata(nanind(length(nanind)-k+1)) = [];
                ydata(nanind(length(nanind)-k+1)) = [];
            end
        end
    end
    allxdata = [allxdata, xdata];
    allydata = [allydata, ydata];
    plot(xdata,ydata,'-o','LineWidth',2)
%     err_std = [this_scbins(i).sdhr_lowsteps, this_scbins(i).sdhr_medwalk, ...
%         this_scbins(i).sdhr_pursteps, this_scbins(i).sdhr_slowwalk, ...
%         this_scbins(i).sdhr_medwalk, this_scbins(i).sdhr_briskwalk, ...
%         this_scbins(i).sdhr_fastamb, this_scbins(i).sdhr_run];
%     errorbar(xdata, ydata, err_std, '-o')

    p = polyfit(xdata,ydata, 1);
    this_scbins(i).slope = p(1);
    this_scbins(i).intercept = p(2);
end

% % Old version: using polyfit on ALL data (no differentiation b/t subjects)
% p = polyfit(allxdata,allydata, 1);
% lin_polyfit_vals(4,1) = p(1);
% lin_polyfit_vals(4,2) = p(2);
% x1 = linspace(0,max(allxdata));
% y1 = polyval(p, x1);
% plot(x1,y1,'Color','k','LineWidth',2)
% TextLocation(['y = ' num2str(p(1)), 'x + ' num2str(p(2))],'Location','southeast');

% New version: taking average slope & intercept to plot LSF line
slope = [];
int = [];
for i = 1:length(this_scbins)
    slope = [slope this_scbins(i).slope];
    int = [int this_scbins(i).intercept];
end
lin_polyfit_vals(4,1) = mean(slope);
lin_polyfit_vals(4,2) = mean(int);
x2 = linspace(0,max(allxdata));
y2 = @(x2) (mean(slope)*x2 + mean(int));
fplot(y2,'-k','LineWidth',2)
TextLocation(['y = ' num2str(mean(slope)),... 
    'x + ' num2str(mean(int))],'Location','southeast');

title('Mean HR at Various SC Thresholds for PAH Patients')
xlabel('SPM')
ylabel('Mean HR (BPM)')
hold off

scbins_PAH = this_scbins;

%% Find which PAH patients are above/below LSF line

% Find which PAH patients are above/below the LSF line.
this_meansc_lowsteps = [];
for i = 1:length(scbins_PAH)
    this_meansc_lowsteps = [this_meansc_lowsteps; scbins_PAH(i).meansc_lowsteps];
end

y_lowsteps = p(1) * mean(this_meansc_lowsteps) + p(2);

indices_below_y_lowsteps = []; % Patients below the LSF line
for i = 1:length(scbins_PAH)
    if scbins_PAH(i).meanhr_lowsteps < y_lowsteps
        indices_below_y_lowsteps = [indices_below_y_lowsteps, i];
    else
        continue
    end
end

for i = 1:length(indices_below_y_lowsteps)
    indices_below_y_lowsteps(i) = PAH_indices(indices_below_y_lowsteps(i));
end

indices_above_y_lowsteps = setdiff(PAH_indices, indices_below_y_lowsteps);


%% Summary of average HR in different SC bins for all subjects

avghealthy_scbins = zeros(1,length(allhealthy_scbins(1,:)));
for i = 1:length(allhealthy_scbins(1,:))
    nanind = find(isnan(allhealthy_scbins(:,i)));
    if isempty(nanind)
        avghealthy_scbins(1,i) = mean(allhealthy_scbins(:,i));
    else
        tempstore = allhealthy_scbins(:,i);
        for j = 1:length(nanind)
            tempstore(nanind(length(nanind)-j+1)) = [];
        end
        avghealthy_scbins(1,i) = mean(tempstore);
    end
end

avgPC_scbins = zeros(1,length(allPC_scbins(1,:)));
for i = 1:length(allPC_scbins(1,:))
    nanind = find(isnan(allPC_scbins(:,i)));
    if isempty(nanind)
        avgPC_scbins(1,i) = mean(allPC_scbins(:,i));
    else
        tempstore = allPC_scbins(:,i);
        for j = 1:length(nanind)
            tempstore(nanind(length(nanind)-j+1)) = [];
        end
        avgPC_scbins(1,i) = mean(tempstore);
    end
end

avgCOPD_scbins = zeros(1,length(allCOPD_scbins(1,:)));
for i = 1:length(allCOPD_scbins(1,:))
    nanind = find(isnan(allCOPD_scbins(:,i)));
    if isempty(nanind)
        avgCOPD_scbins(1,i) = mean(allCOPD_scbins(:,i));
    else
        tempstore = allCOPD_scbins(:,i);
        for j = 1:length(nanind)
            tempstore(nanind(length(nanind)-j+1)) = [];
        end
        avgCOPD_scbins(1,i) = mean(tempstore);
    end
end

avgPAH_scbins = zeros(1,length(allPAH_scbins(1,:)));
for i = 1:length(allPAH_scbins(1,:))
    nanind = find(isnan(allPAH_scbins(:,i)));
    if isempty(nanind)
        avgPAH_scbins(1,i) = mean(allPAH_scbins(:,i));
    else
        tempstore = allPAH_scbins(:,i);
        for j = 1:length(nanind)
            tempstore(nanind(length(nanind)-j+1)) = [];
        end
        avgPAH_scbins(1,i) = mean(tempstore);
    end
end

allavg_scbins = [avghealthy_scbins; avgPC_scbins; avgCOPD_scbins; avgPAH_scbins];

figure(5)
hold on
xdata = [10,30,50,70,90,110,130,150];
plot(xdata, allavg_scbins(1,:),xdata, allavg_scbins(2,:),...
    xdata, allavg_scbins(3,:),xdata, allavg_scbins(4,:),'LineWidth',2)
title('Mean HR at Various SC Thresholds for All Patients')
xlabel('SPM')
ylabel('Mean HR (BPM)')
legend('Healthy', 'PC', 'COPD', 'PAH','Location','best')
hold off

thislegend = ["Healthy Subjects", "PC Patients", "COPD Patients", "PAH Patients"];
figure(6)
hold on
for i = 1:length(lin_polyfit_vals)
    fplot(@(x) lin_polyfit_vals(i,1)*x + lin_polyfit_vals(i,2), [0 150],'LineWidth',2)
end
hold off
title('Mean HR at Various SC Thresholds for All Patients (Linear Fit)')
xlabel('SPM')
ylabel('Mean HR (BPM)')
% TextLocation(['y = ' lin_polyfit_vals(1,1), 'x + ' lin_polyfit_vals(1,2),...
%     ' for ' thislegend(1), 'y = ' lin_polyfit_vals(2,1), 'x + ' lin_polyfit_vals(2,2),...
%     ' for ' thislegend(2), 'y = ' lin_polyfit_vals(3,1), 'x + ' lin_polyfit_vals(3,2),...
%     ' for ' thislegend(3), 'y = ' lin_polyfit_vals(4,1), 'x + ' lin_polyfit_vals(4,2),...
%     ' for ' thislegend(4)],'Location','northwest');
s1 = sprintf('y = (%.2f) x + (%.2f) for Healthy Subjects ',lin_polyfit_vals(1,1),...
    lin_polyfit_vals(1,2));
s2 = sprintf('y = (%.2f) x + (%.2f) for PC Subjects ',lin_polyfit_vals(2,1),...
    lin_polyfit_vals(2,2));
s3 = sprintf('y = (%.2f) x + (%.2f) for COPD Subjects ',lin_polyfit_vals(3,1),...
    lin_polyfit_vals(3,2));
s4 = sprintf('y = (%.2f) x + (%.2f) for PAH Subjects',lin_polyfit_vals(4,1),...
    lin_polyfit_vals(4,2));
TextLocation([s1,newline,s2,newline,s3,newline,s4],'Location','northwest')
legend('Healthy', 'PC', 'COPD', 'PAH','Location','best')

% To plot individual HR vs SC bin
% figure(7)
% hold on
% xdata = [50,70,90,110,130,150];
% group = scbins_PAH; % Input which subject group
% indiv = 7; % Input which row the individual's data is in
% ydata = [group(indiv).meanhr_pursteps, group(indiv).meanhr_slowwalk, ...
%         group(indiv).meanhr_medwalk, group(indiv).meanhr_briskwalk, ...
%         group(indiv).meanhr_fastamb, group(indiv).meanhr_run];
% % plot(xdata, ydata,'-o')
% err_std = [group(indiv).sdhr_pursteps, group(indiv).sdhr_slowwalk, ...
%     group(indiv).sdhr_medwalk, group(indiv).sdhr_briskwalk, ...
%     group(indiv).sdhr_fastamb, group(indiv).sdhr_run];
% errorbar(xdata, ydata, err_std, '-o')
% xlim([40 160])
% title('Mean HR at Various SC Thresholds for PAH09') % Update
% xlabel('SPM')
% ylabel('Mean HR (BPM)')
% hold off

