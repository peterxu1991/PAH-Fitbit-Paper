%% Clear workspace

clear;clc;close all;

%% Organize data for healthy, PC, COPD, and PAH patients

load('Daily_Healthy_PC_110220.mat'); % This structure contains data for healthy, PC, and COPD patients.

% List indices for COPD patients in structure
COPD_indices = [10; 11; 14; 15; 16];

% Create new structures for each patient population.
weekly_raw_COPD = struct([]);

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
% # of SC data points needed in order to count
screq = 10; % Can change this. Minimum is 1 data point.


%% COPD Subjects: calculate avg HR for SC categories

this_weekly_raw = weekly_raw_COPD;

% Create new structures for each patient population.
for i = 1:length(this_weekly_raw)
    
    this_weekly_raw(i).SCbins = struct([]);
    
    for j = 1:length(this_weekly_raw(i).weekly_raw)
        
        SCdata = this_weekly_raw(i).weekly_raw(j).SC;
        HRdata = this_weekly_raw(i).weekly_raw(j).HR;
        this_weekly_raw(i).SCbins(j).Week = this_weekly_raw(i).weekly_raw(j).WeekNo;
        this_weekly_raw(i).SCbins(j).Year = this_weekly_raw(i).weekly_raw(j).YearNo;
        
        for k = 1:length(SCdata)
            
            sc_lowsteps = []; hr_lowsteps = [];
            sc_medsteps = []; hr_medsteps = [];
            sc_pursteps = []; hr_pursteps = [];
            sc_slowwalk = []; hr_slowwalk = [];
            sc_medwalk = []; hr_medwalk = [];
            sc_briskwalk = []; hr_briskwalk = [];
            sc_fastamb = []; hr_fastamb = [];
            sc_run = []; hr_run = [];
            
            lowsteps = find(SCdata > 0 & SCdata < 20);
            medsteps = find(SCdata > 19 & SCdata < 40);
            pursteps = find(SCdata > 39 & SCdata < 60);
            slowwalk = find(SCdata > 59 & SCdata < 80);
            medwalk = find(SCdata > 79 & SCdata < 100);
            briskwalk = find(SCdata > 99 & SCdata < 120);
            fastamb = find(SCdata > 119 & SCdata < 150);
            run = find(SCdata > 149 & SCdata < 250);
            
            if length(lowsteps) > screq - 1
                for l = 1:length(lowsteps)
                    entry = lowsteps(l);
                    sc_lowsteps = [sc_lowsteps; SCdata(entry)];
                    hr_lowsteps = [hr_lowsteps; HRdata(entry)];
                end
                this_weekly_raw(i).SCbins(j).sc_lowsteps = sc_lowsteps;
                this_weekly_raw(i).SCbins(j).meansc_lowsteps = mean(sc_lowsteps);
                this_weekly_raw(i).SCbins(j).hr_lowsteps = hr_lowsteps;
                this_weekly_raw(i).SCbins(j).meanhr_lowsteps = mean(hr_lowsteps);
            else
                this_weekly_raw(i).SCbins(j).meansc_lowsteps = NaN;
                this_weekly_raw(i).SCbins(j).meanhr_lowsteps = NaN;
            end
            
            if length(medsteps) > screq - 1
                for l = 1:length(medsteps)
                    entry = medsteps(l);
                    sc_medsteps = [sc_medsteps;SCdata(entry)];
                    hr_medsteps = [hr_medsteps; HRdata(entry)];
                end
                this_weekly_raw(i).SCbins(j).sc_medsteps = sc_medsteps;
                this_weekly_raw(i).SCbins(j).meansc_medsteps = mean(sc_medsteps);
                this_weekly_raw(i).SCbins(j).hr_medsteps = hr_medsteps;
                this_weekly_raw(i).SCbins(j).meanhr_medsteps = mean(hr_medsteps);
            else
                this_weekly_raw(i).SCbins(j).meansc_medsteps = NaN;
                this_weekly_raw(i).SCbins(j).meanhr_medsteps = NaN;
            end
            
            if length(pursteps) > screq - 1
                for l = 1:length(pursteps)
                    entry = pursteps(l);
                    sc_pursteps = [sc_pursteps;SCdata(entry)];
                    hr_pursteps = [hr_pursteps; HRdata(entry)];
                end
                this_weekly_raw(i).SCbins(j).sc_pursteps = sc_pursteps;
                this_weekly_raw(i).SCbins(j).meansc_pursteps = mean(sc_pursteps);
                this_weekly_raw(i).SCbins(j).hr_pursteps = hr_pursteps;
                this_weekly_raw(i).SCbins(j).meanhr_pursteps = mean(hr_pursteps);
            else
                this_weekly_raw(i).SCbins(j).meansc_pursteps = NaN;
                this_weekly_raw(i).SCbins(j).meanhr_pursteps = NaN;
            end
            
            if length(slowwalk) > screq - 1
                for l = 1:length(slowwalk)
                    entry = slowwalk(l);
                    sc_slowwalk = [sc_slowwalk; SCdata(entry)];
                    hr_slowwalk = [hr_slowwalk; HRdata(entry)];
                end
                this_weekly_raw(i).SCbins(j).sc_slowwalk = sc_slowwalk;
                this_weekly_raw(i).SCbins(j).meansc_slowwalk = mean(sc_slowwalk);
                this_weekly_raw(i).SCbins(j).hr_slowwalk = hr_slowwalk;
                this_weekly_raw(i).SCbins(j).meanhr_slowwalk = mean(hr_slowwalk);
            else
                this_weekly_raw(i).SCbins(j).meansc_slowwalk = NaN;
                this_weekly_raw(i).SCbins(j).meanhr_slowwalk = NaN;
            end
            
            if length(medwalk) > screq - 1
                for l = 1:length(medwalk)
                    entry = medwalk(l);
                    sc_medwalk = [sc_medwalk; SCdata(entry)];
                    hr_medwalk = [hr_medwalk; HRdata(entry)];
                end
                this_weekly_raw(i).SCbins(j).sc_medwalk = sc_medwalk;
                this_weekly_raw(i).SCbins(j).meansc_medwalk = mean(sc_medwalk);
                this_weekly_raw(i).SCbins(j).hr_medwalk = hr_medwalk;
                this_weekly_raw(i).SCbins(j).meanhr_medwalk = mean(hr_medwalk);
            else
                this_weekly_raw(i).SCbins(j).meansc_medwalk = NaN;
                this_weekly_raw(i).SCbins(j).meanhr_medwalk = NaN;
            end
            
            if length(briskwalk) > screq - 1
                for l = 1:length(briskwalk)
                    entry = briskwalk(l);
                    sc_briskwalk = [sc_briskwalk; SCdata(entry)];
                    hr_briskwalk = [hr_briskwalk; HRdata(entry)];
                end
                this_weekly_raw(i).SCbins(j).sc_briskwalk = sc_briskwalk;
                this_weekly_raw(i).SCbins(j).meansc_briskwalk = mean(sc_briskwalk);
                this_weekly_raw(i).SCbins(j).hr_briskwalk = hr_briskwalk;
                this_weekly_raw(i).SCbins(j).meanhr_briskwalk = mean(hr_briskwalk);
            else
                this_weekly_raw(i).SCbins(j).meansc_briskwalk = NaN;
                this_weekly_raw(i).SCbins(j).meanhr_briskwalk = NaN;
            end
            
            if length(fastamb) > screq - 1
                for l = 1:length(fastamb)
                    entry = fastamb(l);
                    sc_fastamb = [sc_fastamb; SCdata(entry)];
                    hr_fastamb = [hr_fastamb; HRdata(entry)];
                end
                this_weekly_raw(i).SCbins(j).sc_fastamb = sc_fastamb;
                this_weekly_raw(i).SCbins(j).meansc_fastamb = mean(sc_fastamb);
                this_weekly_raw(i).SCbins(j).hr_fastamb = hr_fastamb;
                this_weekly_raw(i).SCbins(j).meanhr_fastamb = mean(hr_fastamb);
            else
                this_weekly_raw(i).SCbins(j).meansc_fastamb = NaN;
                this_weekly_raw(i).SCbins(j).meanhr_fastamb = NaN;
            end
            
            if length(run) > screq - 1
                for l = 1:length(run)
                    entry = run(l);
                    sc_run = [sc_run; SCdata(entry)];
                    hr_run = [hr_run; HRdata(entry)];
                end
                this_weekly_raw(i).SCbins(j).sc_run = sc_run;
                this_weekly_raw(i).SCbins(j).meansc_run = mean(sc_run);
                this_weekly_raw(i).SCbins(j).hr_run = hr_run;
                this_weekly_raw(i).SCbins(j).meanhr_run = mean(hr_run);
            else
                this_weekly_raw(i).SCbins(j).meansc_run = NaN;
                this_weekly_raw(i).SCbins(j).meanhr_run = NaN;
            end
            
        end
        
    end
end

weekly_raw_COPD = this_weekly_raw;


%% Slope / intercept using polyfit for COPD Patients

this_weekly_raw = weekly_raw_COPD;
 
for i = 1:length(this_weekly_raw)
    for j = 1:length(this_weekly_raw(i).weekly_raw)
        % Overall slope / intercept
        xdata = [this_weekly_raw(i).SCbins(j).meansc_lowsteps, this_weekly_raw(i).SCbins(j).meansc_medsteps,...
            this_weekly_raw(i).SCbins(j).meansc_pursteps, this_weekly_raw(i).SCbins(j).meansc_slowwalk, ...
            this_weekly_raw(i).SCbins(j).meansc_medwalk, this_weekly_raw(i).SCbins(j).meansc_briskwalk, ...
            this_weekly_raw(i).SCbins(j).meansc_fastamb, this_weekly_raw(i).SCbins(j).meansc_run];
        ydata = [this_weekly_raw(i).SCbins(j).meanhr_lowsteps, this_weekly_raw(i).SCbins(j).meanhr_medsteps,...
            this_weekly_raw(i).SCbins(j).meanhr_pursteps, this_weekly_raw(i).SCbins(j).meanhr_slowwalk, ...
            this_weekly_raw(i).SCbins(j).meanhr_medwalk, this_weekly_raw(i).SCbins(j).meanhr_briskwalk, ...
            this_weekly_raw(i).SCbins(j).meanhr_fastamb, this_weekly_raw(i).SCbins(j).meanhr_run];
        for k = 1:length(xdata)
            nanind = find(isnan(xdata));
            if isempty(nanind)
                continue
            else
                for l = 1:length(nanind)
                    xdata(nanind(length(nanind)-l+1)) = [];
                    ydata(nanind(length(nanind)-l+1)) = [];
                end
            end
        end
        if length(xdata) < 2
            this_weekly_raw(i).SCbins(j).weekly_slope = [];
            this_weekly_raw(i).SCbins(j).weekly_intercept = [];
        else
            p = polyfit(xdata,ydata, 1);
            this_weekly_raw(i).SCbins(j).weekly_slope = p(1);
            this_weekly_raw(i).SCbins(j).weekly_intercept = p(2);
        end
        
        % Region 1 (0 - 60 SPM) slope / intercept
        xdata1 = [this_weekly_raw(i).SCbins(j).meansc_lowsteps, this_weekly_raw(i).SCbins(j).meansc_medsteps,...
            this_weekly_raw(i).SCbins(j).meansc_pursteps];
        ydata1 = [this_weekly_raw(i).SCbins(j).meanhr_lowsteps, this_weekly_raw(i).SCbins(j).meanhr_medsteps,...
            this_weekly_raw(i).SCbins(j).meanhr_pursteps];
        for k = 1:length(xdata1)
            nanind = find(isnan(xdata1));
            if isempty(nanind)
                continue
            else
                for l = 1:length(nanind)
                    xdata1(nanind(length(nanind)-l+1)) = [];
                    ydata1(nanind(length(nanind)-l+1)) = [];
                end
            end
        end
        if length(xdata1) < 2
            this_weekly_raw(i).SCbins(j).weekly_reg1_slope = [];
            this_weekly_raw(i).SCbins(j).weekly_reg1_intercept = [];
        else
            p1 = polyfit(xdata1,ydata1, 1);
            this_weekly_raw(i).SCbins(j).weekly_reg1_slope = p1(1);
            this_weekly_raw(i).SCbins(j).weekly_reg1_intercept = p1(2);
        end
        
        % Region 2 (60+ SPM) slope / intercept
        xdata2 = [this_weekly_raw(i).SCbins(j).meansc_slowwalk, this_weekly_raw(i).SCbins(j).meansc_medwalk,...
            this_weekly_raw(i).SCbins(j).meansc_briskwalk, this_weekly_raw(i).SCbins(j).meansc_fastamb,...
            this_weekly_raw(i).SCbins(j).meansc_run];
        ydata2 = [this_weekly_raw(i).SCbins(j).meanhr_slowwalk, ...
            this_weekly_raw(i).SCbins(j).meanhr_medwalk, this_weekly_raw(i).SCbins(j).meanhr_briskwalk, ...
            this_weekly_raw(i).SCbins(j).meanhr_fastamb, this_weekly_raw(i).SCbins(j).meanhr_run];
        for k = 1:length(xdata2)
            nanind = find(isnan(xdata2));
            if isempty(nanind)
                continue
            else
                for l = 1:length(nanind)
                    xdata2(nanind(length(nanind)-l+1)) = [];
                    ydata2(nanind(length(nanind)-l+1)) = [];
                end
            end
        end
        if length(xdata2) < 2
            this_weekly_raw(i).SCbins(j).weekly_reg2_slope = [];
            this_weekly_raw(i).SCbins(j).weekly_reg2_intercept = [];
        else
            p2 = polyfit(xdata2,ydata2, 1);
            this_weekly_raw(i).SCbins(j).weekly_reg2_slope = p2(1);
            this_weekly_raw(i).SCbins(j).weekly_reg2_intercept = p2(2);
        end 
        
    end
end
 
for i = 1:length(this_weekly_raw)
    this_weekly_raw(i).weekly_slope = []; this_weekly_raw(i).weekly_int = [];
    this_weekly_raw(i).weekly_slope1 = []; this_weekly_raw(i).weekly_int1 = [];
    this_weekly_raw(i).weekly_slope2 = []; this_weekly_raw(i).weekly_int2 = [];
    for j = 1:length(this_weekly_raw(i).SCbins)
        this_weekly_raw(i).weekly_slope = [this_weekly_raw(i).weekly_slope, ...
            this_weekly_raw(i).SCbins(j).weekly_slope];
        this_weekly_raw(i).weekly_int = [this_weekly_raw(i).weekly_int, ...
            this_weekly_raw(i).SCbins(j).weekly_intercept];
        this_weekly_raw(i).weekly_slope1 = [this_weekly_raw(i).weekly_slope1, ...
            this_weekly_raw(i).SCbins(j).weekly_reg1_slope];
        this_weekly_raw(i).weekly_int1 = [this_weekly_raw(i).weekly_int1, ...
            this_weekly_raw(i).SCbins(j).weekly_reg1_intercept];
        this_weekly_raw(i).weekly_slope2 = [this_weekly_raw(i).weekly_slope2, ...
            this_weekly_raw(i).SCbins(j).weekly_reg2_slope];
        this_weekly_raw(i).weekly_int2 = [this_weekly_raw(i).weekly_int2, ...
            this_weekly_raw(i).SCbins(j).weekly_reg2_intercept];
    end
end

for i = 1:length(this_weekly_raw)
    this_weekly_raw(i).meanslope = mean(this_weekly_raw(i).weekly_slope);
    this_weekly_raw(i).meanint = mean(this_weekly_raw(i).weekly_int);
    this_weekly_raw(i).meanslope1 = mean(this_weekly_raw(i).weekly_slope1);
    this_weekly_raw(i).meanint1 = mean(this_weekly_raw(i).weekly_int1);
    this_weekly_raw(i).meanslope2 = mean(this_weekly_raw(i).weekly_slope2);
    this_weekly_raw(i).meanint2 = mean(this_weekly_raw(i).weekly_int2);
end

weekly_raw_COPD = this_weekly_raw;


%% Plot LSF lines for each individual for COPD Patients

figure(1)
hold on
for i = 1:length(COPD_indices)
    y = @(x) (weekly_raw_COPD(i).meanslope*x + weekly_raw_COPD(i).meanint);
    fplot(y,[0,150], 'LineWidth',2)
end
title('LSF Lines for COPD Patients (average weekly LSF values)')
xlabel('SPM')
ylabel('Mean HR (BPM)')
s1 = sprintf('y = (%.2f) x + (%.2f) (COPD 2)',weekly_raw_COPD(1).meanslope,...
    weekly_raw_COPD(1).meanint);
s2 = sprintf('y = (%.2f) x + (%.2f) (COPD 3)',weekly_raw_COPD(2).meanslope,...
    weekly_raw_COPD(2).meanint);
s3 = sprintf('y = (%.2f) x + (%.2f) (COPD 4)',weekly_raw_COPD(3).meanslope,...
    weekly_raw_COPD(3).meanint);
s4 = sprintf('y = (%.2f) x + (%.2f) (COPD 5)',weekly_raw_COPD(4).meanslope,...
    weekly_raw_COPD(4).meanint);
s5 = sprintf('y = (%.2f) x + (%.2f) (COPD 6)',weekly_raw_COPD(5).meanslope,...
    weekly_raw_COPD(5).meanint);
TextLocation([s1,newline,s2,newline,s3,newline,s4,newline,s5],'Location','northwest')
legend('COPD 2','COPD 3','COPD 4','COPD5','COPD6','Location','southeast')
hold off
 

figure(2)
hold on
for i = 1:length(COPD_indices)
    y1 = @(x) (weekly_raw_COPD(i).meanslope1*x + weekly_raw_COPD(i).meanint1);
    fplot(y1,[0,150], 'LineWidth',2)
end
title('LSF Lines for COPD Patients (average weekly LSF values): Region 1')
xlabel('SPM')
ylabel('Mean HR (BPM)')
s1 = sprintf('y = (%.2f) x + (%.2f) (COPD 2)',weekly_raw_COPD(1).meanslope1,...
    weekly_raw_COPD(1).meanint1);
s2 = sprintf('y = (%.2f) x + (%.2f) (COPD 3)',weekly_raw_COPD(2).meanslope1,...
    weekly_raw_COPD(2).meanint1);
s3 = sprintf('y = (%.2f) x + (%.2f) (COPD 4)',weekly_raw_COPD(3).meanslope1,...
    weekly_raw_COPD(3).meanint1);
s4 = sprintf('y = (%.2f) x + (%.2f) (COPD 5)',weekly_raw_COPD(4).meanslope1,...
    weekly_raw_COPD(4).meanint1);
s5 = sprintf('y = (%.2f) x + (%.2f) (COPD 6)',weekly_raw_COPD(5).meanslope1,...
    weekly_raw_COPD(5).meanint1);
TextLocation([s1,newline,s2,newline,s3,newline,s4,newline,s5],'Location','northwest')
legend('COPD 2','COPD 3','COPD 4','COPD5','COPD6','Location','southeast')
hold off
 

figure(3)
hold on
for i = 1:length(COPD_indices)
    y2 = @(x) (weekly_raw_COPD(i).meanslope2*x + weekly_raw_COPD(i).meanint2);
    fplot(y2,[0,150], 'LineWidth',2)
end
title('LSF Lines for COPD Patients (average weekly LSF values): Region 2')
xlabel('SPM')
ylabel('Mean HR (BPM)')
s1 = sprintf('y = (%.2f) x + (%.2f) (COPD 2)',weekly_raw_COPD(1).meanslope2,...
    weekly_raw_COPD(1).meanint2);
s2 = sprintf('y = (%.2f) x + (%.2f) (COPD 3)',weekly_raw_COPD(2).meanslope2,...
    weekly_raw_COPD(2).meanint2);
s3 = sprintf('y = (%.2f) x + (%.2f) (COPD 4)',weekly_raw_COPD(3).meanslope2,...
    weekly_raw_COPD(3).meanint2);
s4 = sprintf('y = (%.2f) x + (%.2f) (COPD 5)',weekly_raw_COPD(4).meanslope2,...
    weekly_raw_COPD(4).meanint2);
s5 = sprintf('y = (%.2f) x + (%.2f) (COPD 6)',weekly_raw_COPD(5).meanslope2,...
    weekly_raw_COPD(5).meanint2);
TextLocation([s1,newline,s2,newline,s3,newline,s4,newline,s5],'Location','southeast')
legend('COPD 2','COPD 3','COPD 4','COPD5','COPD6','Location','north')
hold off


%% Plot average HR in different SC bins (weekly) for COPD Patients

for i = 1:length(weekly_raw_COPD)
    xdata1 = []; ydata1 = [];
    xdata2 = []; ydata2 = [];
    xdata3 = []; ydata3 = [];
    xdata4 = []; ydata4 = [];
    xdata5 = []; ydata5 = [];
    xdata6 = []; ydata6 = [];
    xdata7 = []; ydata7 = [];
    xdata8 = []; ydata8 = [];
    
    for j = 1:length(weekly_raw_COPD(i).SCbins)
        xdata1 = [xdata1; strcat('Week', string(weekly_raw_COPD(i).SCbins(j).Week),...
            ' Year', string(weekly_raw_COPD(i).SCbins(j).Year))];
        ydata1 = [ydata1; weekly_raw_COPD(i).SCbins(j).meanhr_lowsteps];  
    end
    figure(i)
    hold on
    xdata1 = reordercats(categorical(xdata1),xdata1);
    nanind = find(isnan(ydata1));
    if isempty(nanind)
    else
        for k = 1:length(nanind)
            xdata1(nanind(length(nanind)-k+1)) = [];
            ydata1(nanind(length(nanind)-k+1)) = [];
        end
    end
    plot(xdata1,ydata1,'Color',[0.8500, 0.3250, 0.0980],'LineWidth',1)
    
    for j = 1:length(weekly_raw_COPD(i).SCbins)
        xdata2 = [xdata2; strcat('Week', string(weekly_raw_COPD(i).SCbins(j).Week),...
            ' Year', string(weekly_raw_COPD(i).SCbins(j).Year))];
        ydata2 = [ydata2; weekly_raw_COPD(i).SCbins(j).meanhr_medsteps];  
    end
    xdata2 = reordercats(categorical(xdata2),xdata2);
    nanind = find(isnan(ydata2));
    if isempty(nanind)
    else
        for k = 1:length(nanind)
            xdata2(nanind(length(nanind)-k+1)) = [];
            ydata2(nanind(length(nanind)-k+1)) = [];
        end
    end
    plot(xdata2,ydata2,'Color',[0.3010, 0.7450, 0.9330],'LineWidth',1)
    
    for j = 1:length(weekly_raw_COPD(i).SCbins)
        xdata3 = [xdata3; strcat('Week', string(weekly_raw_COPD(i).SCbins(j).Week),...
            ' Year', string(weekly_raw_COPD(i).SCbins(j).Year))];
        ydata3 = [ydata3; weekly_raw_COPD(i).SCbins(j).meanhr_pursteps];  
    end
    xdata3 = reordercats(categorical(xdata3),xdata3);
    nanind = find(isnan(ydata3));
    if isempty(nanind)
    else
        for k = 1:length(nanind)
            xdata3(nanind(length(nanind)-k+1)) = [];
            ydata3(nanind(length(nanind)-k+1)) = [];
        end
    end
    plot(xdata3,ydata3,'Color','m','LineWidth',1)
    
    for j = 1:length(weekly_raw_COPD(i).SCbins)
        xdata4 = [xdata4; strcat('Week', string(weekly_raw_COPD(i).SCbins(j).Week),...
            ' Year', string(weekly_raw_COPD(i).SCbins(j).Year))];
        ydata4 = [ydata4; weekly_raw_COPD(i).SCbins(j).meanhr_slowwalk];  
    end
    xdata4 = reordercats(categorical(xdata4),xdata4);
    nanind = find(isnan(ydata4));
    if isempty(nanind)
    else
        for k = 1:length(nanind)
            xdata4(nanind(length(nanind)-k+1)) = [];
            ydata4(nanind(length(nanind)-k+1)) = [];
        end
    end
    plot(xdata4,ydata4,'Color','c','LineWidth',1)
    
    for j = 1:length(weekly_raw_COPD(i).SCbins)
        xdata5 = [xdata5; strcat('Week', string(weekly_raw_COPD(i).SCbins(j).Week),...
            ' Year', string(weekly_raw_COPD(i).SCbins(j).Year))];
        ydata5 = [ydata5; weekly_raw_COPD(i).SCbins(j).meanhr_medwalk];  
    end
    xdata5 = reordercats(categorical(xdata5),xdata5);
    nanind = find(isnan(ydata5));
    if isempty(nanind)
    else
        for k = 1:length(nanind)
            xdata5(nanind(length(nanind)-k+1)) = [];
            ydata5(nanind(length(nanind)-k+1)) = [];
        end
    end
    plot(xdata5,ydata5,'Color','r','LineWidth',1)
    
    for j = 1:length(weekly_raw_COPD(i).SCbins)
        xdata6 = [xdata6; strcat('Week', string(weekly_raw_COPD(i).SCbins(j).Week),...
            ' Year', string(weekly_raw_COPD(i).SCbins(j).Year))];
        ydata6 = [ydata6; weekly_raw_COPD(i).SCbins(j).meanhr_briskwalk];  
    end
    xdata6 = reordercats(categorical(xdata6),xdata6);
    nanind = find(isnan(ydata6));
    if isempty(nanind)
    else
        for k = 1:length(nanind)
            xdata6(nanind(length(nanind)-k+1)) = [];
            ydata6(nanind(length(nanind)-k+1)) = [];
        end
    end
    plot(xdata6,ydata6,'Color','k','LineWidth',1)
    
    for j = 1:length(weekly_raw_COPD(i).SCbins)
        xdata7 = [xdata7; strcat('Week', string(weekly_raw_COPD(i).SCbins(j).Week),...
            ' Year', string(weekly_raw_COPD(i).SCbins(j).Year))];
        ydata7 = [ydata7; weekly_raw_COPD(i).SCbins(j).meanhr_fastamb];  
    end
    xdata7 = reordercats(categorical(xdata7),xdata7);
    nanind = find(isnan(ydata7));
    if isempty(nanind)
    else
        for k = 1:length(nanind)
            xdata7(nanind(length(nanind)-k+1)) = [];
            ydata7(nanind(length(nanind)-k+1)) = [];
        end
    end
    plot(xdata7,ydata7,'Color','g','LineWidth',1)
    
    for j = 1:length(weekly_raw_COPD(i).SCbins)
        xdata8 = [xdata8; strcat('Week', string(weekly_raw_COPD(i).SCbins(j).Week),...
            ' Year', string(weekly_raw_COPD(i).SCbins(j).Year))];
        ydata8 = [ydata8; weekly_raw_COPD(i).SCbins(j).meanhr_run];  
    end
    xdata8 = reordercats(categorical(xdata8),xdata8);
    nanind = find(isnan(ydata8));
    if isempty(nanind)
    else
        for k = 1:length(nanind)
            xdata8(nanind(length(nanind)-k+1)) = [];
            ydata8(nanind(length(nanind)-k+1)) = [];
        end
    end
    plot(xdata8,ydata8,'Color','b','LineWidth',1)
    
    legend('Low Steps', 'Medium Steps', 'Purposeful Steps', 'Slow Walk', 'Medium Walk', ...
        'Brisk Walk', 'Fast Ambulation', 'Run','Location','eastoutside')
    title(sprintf('Mean HR for %s',string(this_weekly_raw(i).name)))
    xlabel('Time')
    ylabel('Mean HR (BPM)')
    hold off
end
 
 
%% Plot amount of time spent in each SC category

this_weekly_raw = weekly_raw_COPD;
 
% Amount of time (absolute) spent in each SC category (Method 1 - bar graph)
for i = 1:length(this_weekly_raw)
    xdata = []; 
    for j = 1:length(this_weekly_raw(i).SCbins)
        xdata = [xdata; strcat('Week', string(this_weekly_raw(i).SCbins(j).Week),...
                ' Year', string(this_weekly_raw(i).SCbins(j).Year))];
    end
    
    ydata1 = []; ydata2 = []; ydata3 = []; ydata4 = []; 
    ydata5 = []; ydata6 = []; ydata7 = []; ydata8 = [];
    
    exind = [];
    for j = 1:length(this_weekly_raw(i).SCbins)
        if isempty(this_weekly_raw(i).SCbins(j).Week)
            exind = [exind, j];
        else
            continue
        end
    end
    
    if isfield(this_weekly_raw(i).SCbins, 'sc_lowsteps')
        for j = 1:length(this_weekly_raw(i).SCbins)
            ydata1 = [ydata1; length(this_weekly_raw(i).SCbins(j).sc_lowsteps)];
        end
        
        if isempty(exind)
        else
            for k = 1:length(exind)
                ydata1(exind(length(exind)-k+1)) = [];
            end
        end
    end
    
    if isfield(this_weekly_raw(i).SCbins, 'sc_medsteps')
        for j = 1:length(this_weekly_raw(i).SCbins)
            ydata2 = [ydata2; length(this_weekly_raw(i).SCbins(j).sc_medsteps)];
        end
        
        if isempty(exind)
        else
            for k = 1:length(exind)
                ydata2(exind(length(exind)-k+1)) = [];
            end
        end        
    end
    
    if isfield(this_weekly_raw(i).SCbins, 'sc_pursteps')
        for j = 1:length(this_weekly_raw(i).SCbins)
            ydata3 = [ydata3; length(this_weekly_raw(i).SCbins(j).sc_pursteps)];
        end
        
        if isempty(exind)
        else
            for k = 1:length(exind)
                ydata3(exind(length(exind)-k+1)) = [];
            end
        end        
    end
    
    if isfield(this_weekly_raw(i).SCbins, 'sc_slowwalk')
        for j = 1:length(this_weekly_raw(i).SCbins)
            ydata4 = [ydata4; length(this_weekly_raw(i).SCbins(j).sc_slowwalk)];
        end
        
        if isempty(exind)
        else
            for k = 1:length(exind)
                ydata4(exind(length(exind)-k+1)) = [];
            end
        end        
    end
    
    if isfield(this_weekly_raw(i).SCbins, 'sc_medwalk')
        for j = 1:length(this_weekly_raw(i).SCbins)
            ydata5 = [ydata5; length(this_weekly_raw(i).SCbins(j).sc_medwalk)];
        end
        
        if isempty(exind)
        else
            for k = 1:length(exind)
                ydata5(exind(length(exind)-k+1)) = [];
            end
        end        
    end
    
    if isfield(this_weekly_raw(i).SCbins, 'sc_briskwalk')
        for j = 1:length(this_weekly_raw(i).SCbins)
            ydata6 = [ydata6; length(this_weekly_raw(i).SCbins(j).sc_briskwalk)];
        end
        
        if isempty(exind)
        else
            for k = 1:length(exind)
                ydata6(exind(length(exind)-k+1)) = [];
            end
        end        
    end
    
    if isfield(this_weekly_raw(i).SCbins, 'sc_fastamb')
        for j = 1:length(this_weekly_raw(i).SCbins)
            ydata7 = [ydata7; length(this_weekly_raw(i).SCbins(j).sc_fastamb)];
        end
        
        if isempty(exind)
        else
            for k = 1:length(exind)
                ydata7(exind(length(exind)-k+1)) = [];
            end
        end        
    end
    
    if isfield(this_weekly_raw(i).SCbins, 'sc_run')
        for j = 1:length(this_weekly_raw(i).SCbins)
            ydata8 = [ydata8; length(this_weekly_raw(i).SCbins(j).sc_run)];
        end
        
        if isempty(exind)
        else
            for k = 1:length(exind)
                ydata8(exind(length(exind)-k+1)) = [];
            end
        end        
    end
    
    xdata = reordercats(categorical(xdata),xdata);
    ydata = horzcat(ydata1, ydata2, ydata3, ydata4, ydata5, ydata6, ydata7, ydata8);
    
    figure(i)
    bar(xdata,ydata,'stacked')
    legend('Low Steps', 'Medium Steps', 'Purposeful Steps', 'Slow Walk', 'Medium Walk', ...
        'Brisk Walk', 'Fast Ambulation', 'Run','Location','eastoutside')
    title(sprintf('Minutes in Each SC Category for %s',string(this_weekly_raw(i).name)))
    xlabel('Week')
    ylabel('# of minutes in each SC category')
    
end
 
weekly_raw_COPD = this_weekly_raw;
 
 
%% COPD Forecasting: calculate mean, SD, and LSF line

weekthresh = 10; % Minimum 10 weeks data for forecasting analysis.
this_weekly_raw = weekly_raw_COPD;

for i = 1:length(this_weekly_raw)
    if length(this_weekly_raw(i).weekly_raw) < weekthresh
    else
        this_weekly_raw(i).forecast = struct([]);
        numentries = length(daily_Healthy(COPD_indices(i)).hr_mu0) - weekthresh;
        
        % Store mean HR and SD for lowsteps, medsteps, and pursteps.
        meanhr_lowsteps = []; meanhr_medsteps = []; meanhr_pursteps = [];
        for j = 1:length(this_weekly_raw(i).weekly_raw)
            meanhr_lowsteps = [meanhr_lowsteps this_weekly_raw(i).SCbins(j).meanhr_lowsteps];
            meanhr_medsteps = [meanhr_medsteps this_weekly_raw(i).SCbins(j).meanhr_medsteps];
            meanhr_pursteps = [meanhr_pursteps this_weekly_raw(i).SCbins(j).meanhr_pursteps];
        end
        
        for j = 1:numentries
            this_weekly_raw(i).forecast(j).WeekNo = this_weekly_raw(i).weekly_raw(j+weekthresh).WeekNo;
            this_weekly_raw(i).forecast(j).YearNo = this_weekly_raw(i).weekly_raw(j+weekthresh).YearNo;
            
            % Lowsteps
            this_lowsteps = meanhr_lowsteps(1:(weekthresh+j-1));
            nanind = find(isnan(this_lowsteps));
            if isempty(nanind)
            else
                for l = 1:length(nanind)
                    this_lowsteps(nanind(length(nanind)-l+1)) = [];
                end
            end
            this_weekly_raw(i).forecast(j).moving_mean_lowsteps = mean(this_lowsteps);
            this_weekly_raw(i).forecast(j).moving_SD_lowsteps = std(this_lowsteps);
            
            % Medsteps
            this_medsteps = meanhr_medsteps(1:(weekthresh+j-1));
            nanind = find(isnan(this_medsteps));
            if isempty(nanind)
            else
                for l = 1:length(nanind)
                    this_medsteps(nanind(length(nanind)-l+1)) = [];
                end
            end
            this_weekly_raw(i).forecast(j).moving_mean_medsteps = mean(this_medsteps);
            this_weekly_raw(i).forecast(j).moving_SD_medsteps = std(this_medsteps);
            
            % Pursteps
            this_pursteps = meanhr_pursteps(1:(weekthresh+j-1));
            nanind = find(isnan(this_pursteps));
            if isempty(nanind)
            else
                for l = 1:length(nanind)
                    this_pursteps(nanind(length(nanind)-l+1)) = [];
                end
            end
            this_weekly_raw(i).forecast(j).moving_mean_pursteps = mean(this_pursteps);
            this_weekly_raw(i).forecast(j).moving_SD_pursteps = std(this_pursteps);
        end
        
        % Store LSF mean and intercept for lowsteps, medsteps, and pursteps.
        for j = 1:numentries
            
            % Lowsteps
            xdata = 1:1:weekthresh+j-1;
            ydata1 = meanhr_lowsteps(1:(weekthresh+j-1));
            nanind = find(isnan(ydata1));
            if isempty(nanind)
            else
                for l = 1:length(nanind)
                    xdata(nanind(length(nanind)-l+1)) = [];
                    ydata1(nanind(length(nanind)-l+1)) = [];
                end
            end
            p = polyfit(xdata, ydata1, 1);
            this_weekly_raw(i).forecast(j).moving_slope_lowsteps = p(1);
            this_weekly_raw(i).forecast(j).moving_int_lowsteps = p(2);
            
            % Medsteps
            xdata = 1:1:weekthresh+j-1;
            ydata2 = meanhr_medsteps(1:(weekthresh+j-1));
            nanind = find(isnan(ydata2));
            if isempty(nanind)
            else
                for l = 1:length(nanind)
                    xdata(nanind(length(nanind)-l+1)) = [];
                    ydata2(nanind(length(nanind)-l+1)) = [];
                end
            end
            p = polyfit(xdata, ydata2, 1);
            this_weekly_raw(i).forecast(j).moving_slope_medsteps = p(1);
            this_weekly_raw(i).forecast(j).moving_int_medsteps = p(2);
            
            % Pursteps
            xdata = 1:1:weekthresh+j-1;
            ydata3 = meanhr_pursteps(1:(weekthresh+j-1));
            nanind = find(isnan(ydata3));
            if isempty(nanind)
            else
                for l = 1:length(nanind)
                    xdata(nanind(length(nanind)-l+1)) = [];
                    ydata3(nanind(length(nanind)-l+1)) = [];
                end
            end
            p = polyfit(xdata, ydata3, 1);
            this_weekly_raw(i).forecast(j).moving_slope_pursteps = p(1);
            this_weekly_raw(i).forecast(j).moving_int_pursteps = p(2);
        end
    end
end
 
weekly_raw_COPD = this_weekly_raw;
 
 
%% Plot Forecast predictions using mean & SD of lowest 3 SC categories

% How many standard deviations above/below the mean HR is acceptable?
n = 2;
this_weekly_raw = weekly_raw_COPD;
 
% Plot 1st lowest SC category (lowsteps)
for i = 1:length(this_weekly_raw)
    if isempty(this_weekly_raw(i).forecast)
    else
        xdata1 = []; ydata1 = [];
        xdata1p = []; ydata1p = [];ydata1p_sd = [];

        for j = 1:length(this_weekly_raw(i).SCbins)
            xdata1 = [xdata1; strcat('Week', string(this_weekly_raw(i).SCbins(j).Week),...
                ' Year', string(this_weekly_raw(i).SCbins(j).Year))];
            ydata1 = [ydata1; this_weekly_raw(i).SCbins(j).meanhr_lowsteps];  
        end
        figure(i)
        hold on
        xdata1 = reordercats(categorical(xdata1),xdata1);
        plot(xdata1,ydata1,'o-','Color',[0.8500, 0.3250, 0.0980],'LineWidth',1)

        % Predicted HR values & standard deviation bands
        for j = 1:length(this_weekly_raw(i).forecast)
            xdata1p = [xdata1p; strcat('Week', string(this_weekly_raw(i).forecast(j).WeekNo),...
                ' Year', string(this_weekly_raw(i).forecast(j).YearNo))];
            ydata1p = [ydata1p; this_weekly_raw(i).forecast(j).moving_mean_lowsteps];
            ydata1p_sd = [ydata1p_sd; this_weekly_raw(i).forecast(j).moving_SD_lowsteps];
        end
        xdata1p = reordercats(categorical(xdata1p),xdata1p);
        meansd = mean(ydata1p_sd);

        plot(xdata1p,ydata1p,'o--','Color','k','LineWidth',1) 
        plot(xdata1p,ydata1p + n*ydata1p_sd,'--','Color','k','LineWidth',1)
        plot(xdata1p,ydata1p - n*ydata1p_sd,'--','Color','k','LineWidth',1)

        TextLocation(['Average distance of', newline, 'upper/lower limit', newline,... 
            'from predicted HR', newline, 'line = ' num2str(sprintf('%.3g',n*meansd))],...
            'Location','northeastoutside');

    %     % Failed attempts to shade region between SD lines
    %     x2 = [xdata1p, fliplr(xdata1p)];
    %     inbetween = [ydata1p - ydata1p_sd, fliplr(ydata1p + ydata1p_sd)];
    %     fill(x2, inbetween, 'g');
    %     patch(x2, inbetween, 'g')
    %     stdshade(ydata1p_sd,0.5,'g');

        legend('Low Steps','Predicted Low Steps','Location','eastoutside')
        title(sprintf('Mean HR for %s',string(this_weekly_raw(i).name)))
        xlabel('Time')
        ylabel('Mean HR (BPM)')
        hold off
    end
end
 
% Plot 2nd lowest SC category (medsteps)
for i = 1:length(this_weekly_raw)
    if isempty(this_weekly_raw(i).forecast)
    else
        xdata2 = []; ydata2 = [];
        xdata2p = []; ydata2p = [];ydata2p_sd = [];

        for j = 1:length(this_weekly_raw(i).SCbins)
            xdata2 = [xdata2; strcat('Week', string(this_weekly_raw(i).SCbins(j).Week),...
                ' Year', string(this_weekly_raw(i).SCbins(j).Year))];
            ydata2 = [ydata2; this_weekly_raw(i).SCbins(j).meanhr_medsteps];  
        end
        figure(i+length(this_weekly_raw))
        hold on
        xdata2 = reordercats(categorical(xdata2),xdata2);
        plot(xdata2,ydata2,'o-','Color','b','LineWidth',1)

        % Predicted HR values & standard deviation bands
        for j = 1:length(this_weekly_raw(i).forecast)
            xdata2p = [xdata2p; strcat('Week', string(this_weekly_raw(i).forecast(j).WeekNo),...
                ' Year', string(this_weekly_raw(i).forecast(j).YearNo))];
            ydata2p = [ydata2p; this_weekly_raw(i).forecast(j).moving_mean_medsteps];
            ydata2p_sd = [ydata2p_sd; this_weekly_raw(i).forecast(j).moving_SD_medsteps];
        end
        xdata2p = reordercats(categorical(xdata2p),xdata2p);
        meansd = mean(ydata2p_sd);

        plot(xdata2p,ydata2p,'o--','Color','k','LineWidth',1) 
        plot(xdata2p,ydata2p + n*ydata2p_sd,'--','Color','k','LineWidth',1)
        plot(xdata2p,ydata2p - n*ydata2p_sd,'--','Color','k','LineWidth',1)

        TextLocation(['Average distance of', newline, 'upper/lower limit', newline,... 
            'from predicted HR', newline, 'line = ' num2str(sprintf('%.3g',n*meansd))],...
            'Location','northeastoutside');

        legend('Medium (Sporadic) Steps','Predicted Medium Steps','Location','eastoutside')
        title(sprintf('Mean HR for %s',string(this_weekly_raw(i).name)))
        xlabel('Time')
        ylabel('Mean HR (BPM)')
        hold off
    end   
end
 
% Plot 3rd lowest SC category (pursteps)
for i = 1:length(this_weekly_raw)
    if isempty(this_weekly_raw(i).forecast)
    else
        xdata3 = []; ydata3 = [];
        xdata3p = []; ydata3p = [];ydata3p_sd = [];

        for j = 1:length(this_weekly_raw(i).SCbins)
            xdata3 = [xdata3; strcat('Week', string(this_weekly_raw(i).SCbins(j).Week),...
                ' Year', string(this_weekly_raw(i).SCbins(j).Year))];
            ydata3 = [ydata3; this_weekly_raw(i).SCbins(j).meanhr_pursteps];  
        end
        figure(i+2*length(this_weekly_raw))
        hold on
        xdata3 = reordercats(categorical(xdata3),xdata3);
        plot(xdata3,ydata3,'o-','Color','m','LineWidth',1)

        % Predicted HR values & standard deviation bands
        for j = 1:length(this_weekly_raw(i).forecast)
            xdata3p = [xdata3p; strcat('Week', string(this_weekly_raw(i).forecast(j).WeekNo),...
                ' Year', string(this_weekly_raw(i).forecast(j).YearNo))];
            ydata3p = [ydata3p; this_weekly_raw(i).forecast(j).moving_mean_pursteps];
            ydata3p_sd = [ydata3p_sd; this_weekly_raw(i).forecast(j).moving_SD_pursteps];
        end
        xdata3p = reordercats(categorical(xdata3p),xdata3p);
        meansd = mean(ydata3p_sd);

        plot(xdata3p,ydata3p,'o--','Color','k','LineWidth',1) 
        plot(xdata3p,ydata3p + n*ydata3p_sd,'--','Color','k','LineWidth',1)
        plot(xdata3p,ydata3p - n*ydata3p_sd,'--','Color','k','LineWidth',1)

        TextLocation(['Average distance of', newline, 'upper/lower limit', newline,... 
            'from predicted HR', newline, 'line = ' num2str(sprintf('%.3g',n*meansd))],...
            'Location','northeastoutside');

        legend('Purposeful Steps','Predicted Purposeful Steps','Location','eastoutside')
        title(sprintf('Mean HR for %s',string(this_weekly_raw(i).name)))
        xlabel('Time')
        ylabel('Mean HR (BPM)')
        hold off
    end 
end

weekly_raw_COPD = this_weekly_raw;

 
%% Plot Forecast predictions using LSF line of lowest 3 SC categories

this_weekly_raw = weekly_raw_COPD;
 
% Plot 1st lowest SC category (lowsteps)
for i = 1:length(this_weekly_raw)
    if isempty(this_weekly_raw(i).forecast)
    else
        xdata1 = []; ydata1 = [];
        xdata1p = []; ydata1p = [];ydata1p_sd = [];
 
        for j = 1:length(this_weekly_raw(i).SCbins)
            xdata1 = [xdata1; strcat('Week', string(this_weekly_raw(i).SCbins(j).Week),...
                ' Year', string(this_weekly_raw(i).SCbins(j).Year))];
            ydata1 = [ydata1; this_weekly_raw(i).SCbins(j).meanhr_lowsteps];  
        end
        figure(i)
        hold on
        xdata1 = reordercats(categorical(xdata1),xdata1);
        plot(xdata1,ydata1,'o-','Color',[0.8500, 0.3250, 0.0980],'LineWidth',1)
 
        % LSF line
        x = linspace(0,length(xdata1));
        lastent = length(this_weekly_raw(i).forecast);
        slope = this_weekly_raw(i).forecast(lastent).moving_slope_lowsteps;
        int = this_weekly_raw(i).forecast(lastent).moving_int_lowsteps;
        y = @(x) (slope*x + int);
        fplot(y,'-k','LineWidth',2)
        TextLocation(['y = ' num2str(slope),... 
            'x + ' num2str(int)],'Location','northeastoutside');
 
        legend('Low Steps','LSF Line for Low Steps','Location','eastoutside')
        title(sprintf('Mean HR for %s',string(this_weekly_raw(i).name)))
        xlabel('Time')
        ylabel('Mean HR (BPM)')
        hold off
    end
end
 
% Plot 2nd lowest SC category (medsteps)
for i = 1:length(this_weekly_raw)
    if isempty(this_weekly_raw(i).forecast)
    else
        xdata2 = []; ydata2 = [];
        xdata2p = []; ydata2p = [];ydata2p_sd = [];
 
        for j = 1:length(this_weekly_raw(i).SCbins)
            xdata2 = [xdata2; strcat('Week', string(this_weekly_raw(i).SCbins(j).Week),...
                ' Year', string(this_weekly_raw(i).SCbins(j).Year))];
            ydata2 = [ydata2; this_weekly_raw(i).SCbins(j).meanhr_medsteps];  
        end
        figure(i+length(this_weekly_raw))
        hold on
        xdata2 = reordercats(categorical(xdata2),xdata2);
        plot(xdata2,ydata2,'o-','Color','b','LineWidth',1)
 
        % LSF line
        x = linspace(0,length(xdata2));
        lastent = length(this_weekly_raw(i).forecast);
        slope = this_weekly_raw(i).forecast(lastent).moving_slope_medsteps;
        int = this_weekly_raw(i).forecast(lastent).moving_int_medsteps;
        y = @(x) (slope*x + int);
        fplot(y,'-k','LineWidth',2)
        TextLocation(['y = ' num2str(slope),... 
            'x + ' num2str(int)],'Location','northeastoutside');
 
        legend('Medium (Sporadic) Steps','LSF Line for Medium Steps','Location','eastoutside')
        title(sprintf('Mean HR for %s',string(this_weekly_raw(i).name)))
        xlabel('Time')
        ylabel('Mean HR (BPM)')
        hold off
    end
end

% Plot 3rd lowest SC category (pursteps)
for i = 1:length(this_weekly_raw)
    if isempty(this_weekly_raw(i).forecast)
    else
        xdata3 = []; ydata3 = [];
        xdata3p = []; ydata3p = [];ydata3p_sd = [];
 
        for j = 1:length(this_weekly_raw(i).SCbins)
            xdata3 = [xdata3; strcat('Week', string(this_weekly_raw(i).SCbins(j).Week),...
                ' Year', string(this_weekly_raw(i).SCbins(j).Year))];
            ydata3 = [ydata3; this_weekly_raw(i).SCbins(j).meanhr_pursteps];  
        end
        figure(i+2*length(this_weekly_raw))
        hold on
        xdata3 = reordercats(categorical(xdata3),xdata3);
        plot(xdata3,ydata3,'o-','Color','m','LineWidth',1)
 
        % LSF line
        x = linspace(0,length(xdata3));
        lastent = length(this_weekly_raw(i).forecast);
        slope = this_weekly_raw(i).forecast(lastent).moving_slope_pursteps;
        int = this_weekly_raw(i).forecast(lastent).moving_int_pursteps;
        y = @(x) (slope*x + int);
        fplot(y,'-k','LineWidth',2)
        TextLocation(['y = ' num2str(slope),... 
            'x + ' num2str(int)],'Location','northeastoutside');
 
        legend('Purposeful Steps','LSF Line for Purposeful Steps','Location','eastoutside')
        title(sprintf('Mean HR for %s',string(this_weekly_raw(i).name)))
        xlabel('Time')
        ylabel('Mean HR (BPM)')
        hold off
    end
end

weekly_raw_COPD = this_weekly_raw;

