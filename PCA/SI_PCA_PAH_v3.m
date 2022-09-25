%%% Use this code for generating PCA plots for only PAH patients
%%% This code takes X # of iterations of analysis before generating plots.

%%% Statistics of Fitbit step & Heart rate

clear;clc;close all;

%% read data
load('daily_PAH_082520.mat');
load('PAH_PCA_sustained param.mat');
load('PAH_PCA_sedentary.mat');
PAH_clindata = readmatrix('PAH_Clin Data.xlsx');
PAH_ratios = readmatrix('AK_PAH_Ratios.csv');
PAH_compliance_usage = readmatrix('AK_PAH_compliance_usage.xlsx');
PAH_trunc_clindata = PAH_clindata(1:30, [1,82:85]);

% List indices for PAH patients in structure
PAH_ex_indices = [6; 16; 18; 24; 25; 29]; % Patient indices to exclude
PAH_indices = 1:1:30; % Array with all PAH patients (even patients to exclude)

% List of patient indices to be included in analysis
for i = 1:length(PAH_ex_indices)
    PAH_indices(PAH_indices == PAH_ex_indices(i)) = [];
end

%%% General Parameters %%%
thresh = 10; % Threshold for number of weeks in PCA per participant.
iterations = 10; % # of iterations of PCA to perform before plotting.

%% list the variables needed and convert the overall structure to cell
% without usage
var_name = {'hr_mu0','hr_sig0','hr_sk0','hr_ks0','hr_sig1','hr_sk1','hr_ks1',...
    'sc_sk','sc_ks','slope','area','fA1','fA3',...
    'mean_top_sc','mean_top_hr','MWT_sc','MWT_hr',...
    'facti','day_pct','night_pct','compliance','missing_hours',...
    'end_expdist_mu','freq_sustained','age','weight','height',...
    'compliance2','maxofftime','ratio6MWT'};

var_name_trunc1 = {'hr_mu0','hr_sig0','hr_sk0','hr_sig1','hr_sk1',...
    'sc_mu','sc_sig','sc_sk'};

var_name_trunc1a = {'hr_mu0','hr_sk0','hr_sig1','hr_sk1',...
    'sc_mu','sc_sig','sc_sk'};

var_name_trunc1b = {'hr_mu0','hr_sk0','hr_sig1','hr_sk1',...
    'sc_mu','sc_sig','sc_sk','sedentary'};

var_name_trunc1c = {'hr_mu0','hr_sk0','hr_sig1','hr_sk1',...
    'sc_mu','sc_sig','sedentary'};

var_name_trunc1d = {'hr_mu0','hr_sig1','sc_mu','sc_sig','sedentary'};

var_name_trunc1e = {'hr_mu0','hr_sig1','hr_sk1',...
    'sc_mu','sc_sig','sedentary'};

var_name_trunc2 = {'hr_mu0','hr_sig0','hr_sk0','hr_sig1','hr_sk1',...
    'sc_mu','sc_sig','sc_sk','slope','area','facti','age','weight','height'};

var_name_derived = {'slope','area','fA1','fA3','MWT_sc','MWT_hr',...
    'facti','end_expdist_mu','freq_sustained','age','weight','height',...
    'compliance2','maxofftime','ratio6MWT'};

% Match up variables with corresponding indices
hr0_ind = [3,4,5,7];
hr1_ind = [8,9,10,12]; % heart rate data distribution: mu & sigma & skewness & Kolmogorov Smirnov test statistics (ks)
sc_ind = [13,14,15,17]; % step count data distribution: mu & sigma & skewness & Kolmogorov Smirnov test statistics (ks)
hr_sc_ind = [18,19]; % mapping statistics: slope & intercept&area
area_ind = 20; % area
act_pct_ind = 22; % fraction of total time at each activity level
wear_pct_ind = [23,24,25,26,27];%fraction of active time (SC>0),day time&nighttime percentage,compliance,missing hours
hour_step_ind = [28,29,30];
max_sc_hr_ind = [31,32,33,34];
peak_pf_ind = [35,36];
MWT_ind = [37,38];
Pred_SC_HR_ind = [39,40];
sustained_walk_ind = [41,42,43];
sedentary_ind = [44];
demographics_ind = [1,45,46];
usage_ratio_ind = [47,48,49];

% Store all PAH data in daily_PAH structure.
new_daily_PAH = struct([]);
for i = 1:length(PAH_indices)
    new_daily_PAH = [new_daily_PAH; daily_PAH(PAH_indices(i))];
end

% Add sustained walking parameters & weight & height
for i = 1:length(PAH_indices)
    new_daily_PAH(i).end_expdist_mu = PAH_pca_struct(i).end_expdist_mu;
    new_daily_PAH(i).freq_sustained = PAH_pca_struct(i).freq;
    new_daily_PAH(i).int_expdist_mu = PAH_pca_struct(i).int_expdist_mu;
    new_daily_PAH(i).weekly_mean0SCpct = PCA_sedentary(i).weekly_mean0SCpct;
    new_daily_PAH(i).weight = PAH_trunc_clindata(i,2);
    new_daily_PAH(i).height = sqrt(PAH_trunc_clindata(i,2)/PAH_trunc_clindata(i,4));
end

% Add compliance and maximum off-time values
for i = 1:length(PAH_indices)
    temp1 = []; temp2 = [];
    for j = 1:length(PAH_compliance_usage)
        if PAH_compliance_usage(j,1) == PAH_indices(i)
            temp1 = [temp1; PAH_compliance_usage(j,2)];
            temp2 = [temp2; PAH_compliance_usage(j,3)];
        end
    end
    new_daily_PAH(i).AKcompliance = temp1;
    new_daily_PAH(i).maxofftime = temp2;
end

% Add 6MWT ratios (clarify with AK what this is)
for i = 1:length(PAH_indices)
    temp = [];
    for j = 1:size(PAH_ratios,1)
        if PAH_ratios(j,1) == PAH_indices(i)
            temp = PAH_ratios(j,2:end);
        end
    end
    temp = temp(~isnan(temp));
    new_daily_PAH(i).ratio6MWT = temp';
end

% Find PAH participants who have at least {thresh} weeks of data.
PAH_include_indices = [];
for i = 1:length(new_daily_PAH)
    if length(new_daily_PAH(i).hr_mu0) > thresh - 1
        PAH_include_indices = [PAH_include_indices i];
    end
end

% Store PAH participants who have at least {thresh} weeks of data in new structure.
daily_PAH_include = struct([]);
for i = 1:length(PAH_include_indices)
    daily_PAH_include = [daily_PAH_include; new_daily_PAH(PAH_include_indices(i))];
end


%% Perform X iterations of PCA

% sum_score_all1 = zeros(length(daily_PAH_include)*thresh, iterations); 
% sum_score_all2 = zeros(length(daily_PAH_include)*thresh, iterations);
% sum_explained = zeros(length(var_name), iterations);
% sum_coeff_all1 = zeros(length(var_name), iterations); 
% sum_coeff_all2 = zeros(length(var_name), iterations);

numtrial = 100; % # of trials to run; change as needed.

variance_captured = zeros(numtrial,1);

for k = 1:numtrial
    rng(k); % 15 is the one we are using in SI_PCA_PAH_v2.m 
    % Structure with randomly selected {thresh} weeks for selected PAH participants
    daily_PAH_include_trun = struct([]);
    
    for i = 1:length(daily_PAH_include)
        thislength = length(daily_PAH_include(i).AKcompliance); % FIXME: using AKcompliance bc some have less weeks...
        selectweeks = randperm(thislength, thresh);
        daily_PAH_include_trun(i).age = daily_PAH_include(i).age;
        daily_PAH_include_trun(i).name = daily_PAH_include(i).name;
        daily_PAH_include_trun(i).hr_mu0 = [];
        daily_PAH_include_trun(i).hr_sig0 = [];
        daily_PAH_include_trun(i).hr_sk0 = [];
        daily_PAH_include_trun(i).hr_ku0 = [];
        daily_PAH_include_trun(i).hr_ks0 = [];
        daily_PAH_include_trun(i).hr_mu1 = [];
        daily_PAH_include_trun(i).hr_sig1 = [];
        daily_PAH_include_trun(i).hr_sk1 = [];
        daily_PAH_include_trun(i).hr_ku1 = [];
        daily_PAH_include_trun(i).hr_ks1 = [];
        daily_PAH_include_trun(i).sc_mu = [];
        daily_PAH_include_trun(i).sc_sig = [];
        daily_PAH_include_trun(i).sc_sk = [];
        daily_PAH_include_trun(i).sc_ku = [];
        daily_PAH_include_trun(i).sc_ks = [];
        daily_PAH_include_trun(i).slope = [];
        daily_PAH_include_trun(i).intercept = [];
        daily_PAH_include_trun(i).area = [];
        daily_PAH_include_trun(i).act_abs = [];
        daily_PAH_include_trun(i).act_pct = [];
        daily_PAH_include_trun(i).facti = [];
        daily_PAH_include_trun(i).day_pct = [];
        daily_PAH_include_trun(i).night_pct = [];
        daily_PAH_include_trun(i).compliance = [];
        daily_PAH_include_trun(i).missing_hours = [];
        daily_PAH_include_trun(i).area_step_week = [];
        daily_PAH_include_trun(i).max_hour_step = [];
        daily_PAH_include_trun(i).core_hour1 = [];
        daily_PAH_include_trun(i).max_step = [];
        daily_PAH_include_trun(i).hr_max_sc = [];
        daily_PAH_include_trun(i).max_hr = [];
        daily_PAH_include_trun(i).sc_max_hr = [];
        daily_PAH_include_trun(i).mean_top_sc = [];
        daily_PAH_include_trun(i).mean_top_hr = [];
        daily_PAH_include_trun(i).MWT_sc = [];
        daily_PAH_include_trun(i).MWT_hr = [];
        daily_PAH_include_trun(i).SC170 = [];
        daily_PAH_include_trun(i).HR170 = [];
        daily_PAH_include_trun(i).end_expdist_mu = [];
        daily_PAH_include_trun(i).freq_sustained = [];
        daily_PAH_include_trun(i).int_expdist_mu = [];
        daily_PAH_include_trun(i).weekly_mean0SCpct = [];
        daily_PAH_include_trun(i).weight = daily_PAH_include(i).weight;
        daily_PAH_include_trun(i).height = daily_PAH_include(i).height;
        daily_PAH_include_trun(i).AKcompliance = [];
        daily_PAH_include_trun(i).maxofftime = [];
        daily_PAH_include_trun(i).ratio6MWT = [];
        % daily_PAH_include_trun(i).weekly_raw = struct([]);
        for j = 1:length(selectweeks)
            daily_PAH_include_trun(i).hr_mu0 = [daily_PAH_include_trun(i).hr_mu0; daily_PAH_include(i).hr_mu0(selectweeks(j))];
            daily_PAH_include_trun(i).hr_sig0 = [daily_PAH_include_trun(i).hr_sig0; daily_PAH_include(i).hr_sig0(selectweeks(j))];
            daily_PAH_include_trun(i).hr_sk0 = [daily_PAH_include_trun(i).hr_sk0; daily_PAH_include(i).hr_sk0(selectweeks(j))];
            daily_PAH_include_trun(i).hr_ku0 = [daily_PAH_include_trun(i).hr_ku0; daily_PAH_include(i).hr_ku0(selectweeks(j))];
            daily_PAH_include_trun(i).hr_ks0 = [daily_PAH_include_trun(i).hr_ks0; daily_PAH_include(i).hr_ks0(selectweeks(j))];
            daily_PAH_include_trun(i).hr_mu1 = [daily_PAH_include_trun(i).hr_mu1; daily_PAH_include(i).hr_mu1(selectweeks(j))];
            daily_PAH_include_trun(i).hr_sig1 = [daily_PAH_include_trun(i).hr_sig1; daily_PAH_include(i).hr_sig1(selectweeks(j))];
            daily_PAH_include_trun(i).hr_sk1 = [daily_PAH_include_trun(i).hr_sk1; daily_PAH_include(i).hr_sk1(selectweeks(j))];
            daily_PAH_include_trun(i).hr_ku1 = [daily_PAH_include_trun(i).hr_ku1; daily_PAH_include(i).hr_ku1(selectweeks(j))];
            daily_PAH_include_trun(i).hr_ks1 = [daily_PAH_include_trun(i).hr_ks1; daily_PAH_include(i).hr_ks1(selectweeks(j))];
            daily_PAH_include_trun(i).sc_mu = [daily_PAH_include_trun(i).sc_mu; daily_PAH_include(i).sc_mu(selectweeks(j))];
            daily_PAH_include_trun(i).sc_sig = [daily_PAH_include_trun(i).sc_sig; daily_PAH_include(i).sc_sig(selectweeks(j))];
            daily_PAH_include_trun(i).sc_sk = [daily_PAH_include_trun(i).sc_sk; daily_PAH_include(i).sc_sk(selectweeks(j))];
            daily_PAH_include_trun(i).sc_ku = [daily_PAH_include_trun(i).sc_ku; daily_PAH_include(i).sc_ku(selectweeks(j))];
            daily_PAH_include_trun(i).sc_ks = [daily_PAH_include_trun(i).sc_ks; daily_PAH_include(i).sc_ks(selectweeks(j))];
            daily_PAH_include_trun(i).slope = [daily_PAH_include_trun(i).slope; daily_PAH_include(i).slope(selectweeks(j))];
            daily_PAH_include_trun(i).intercept = [daily_PAH_include_trun(i).intercept; daily_PAH_include(i).intercept(selectweeks(j))];
            daily_PAH_include_trun(i).area = [daily_PAH_include_trun(i).area; daily_PAH_include(i).area(selectweeks(j),:)];
            daily_PAH_include_trun(i).act_abs = [daily_PAH_include_trun(i).act_abs; daily_PAH_include(i).act_abs(selectweeks(j),:)];
            daily_PAH_include_trun(i).act_pct = [daily_PAH_include_trun(i).act_pct; daily_PAH_include(i).act_pct(selectweeks(j),:)];
            daily_PAH_include_trun(i).facti = [daily_PAH_include_trun(i).facti; daily_PAH_include(i).facti(selectweeks(j))];
            daily_PAH_include_trun(i).day_pct = [daily_PAH_include_trun(i).day_pct; daily_PAH_include(i).day_pct(selectweeks(j))];
            daily_PAH_include_trun(i).night_pct = [daily_PAH_include_trun(i).night_pct; daily_PAH_include(i).night_pct(selectweeks(j))];
            daily_PAH_include_trun(i).compliance = [daily_PAH_include_trun(i).compliance; daily_PAH_include(i).compliance(selectweeks(j))];
            daily_PAH_include_trun(i).missing_hours = [daily_PAH_include_trun(i).missing_hours; daily_PAH_include(i).missing_hours(selectweeks(j))];
            daily_PAH_include_trun(i).area_step_week = [daily_PAH_include_trun(i).area_step_week; daily_PAH_include(i).area_step_week(selectweeks(j))];
            daily_PAH_include_trun(i).max_hour_step = [daily_PAH_include_trun(i).max_hour_step; daily_PAH_include(i).max_hour_step(selectweeks(j))];
            daily_PAH_include_trun(i).core_hour1 = [daily_PAH_include_trun(i).core_hour1; daily_PAH_include(i).core_hour1(selectweeks(j))];
            daily_PAH_include_trun(i).max_step = [daily_PAH_include_trun(i).max_step; daily_PAH_include(i).max_step(selectweeks(j))];
            daily_PAH_include_trun(i).hr_max_sc = [daily_PAH_include_trun(i).hr_max_sc; daily_PAH_include(i).hr_max_sc(selectweeks(j))];
            daily_PAH_include_trun(i).max_hr = [daily_PAH_include_trun(i).max_hr; daily_PAH_include(i).max_hr(selectweeks(j))];
            daily_PAH_include_trun(i).sc_max_hr = [daily_PAH_include_trun(i).sc_max_hr; daily_PAH_include(i).sc_max_hr(selectweeks(j))];
            daily_PAH_include_trun(i).mean_top_sc = [daily_PAH_include_trun(i).mean_top_sc; daily_PAH_include(i).mean_top_sc(selectweeks(j))];
            daily_PAH_include_trun(i).mean_top_hr = [daily_PAH_include_trun(i).mean_top_hr; daily_PAH_include(i).mean_top_hr(selectweeks(j))];
            daily_PAH_include_trun(i).MWT_sc = [daily_PAH_include_trun(i).MWT_sc; daily_PAH_include(i).MWT_sc(selectweeks(j))];
            daily_PAH_include_trun(i).MWT_hr = [daily_PAH_include_trun(i).MWT_hr; daily_PAH_include(i).MWT_hr(selectweeks(j))];
            daily_PAH_include_trun(i).SC170 = [daily_PAH_include_trun(i).SC170; daily_PAH_include(i).SC170(selectweeks(j))];
            daily_PAH_include_trun(i).HR170 = [daily_PAH_include_trun(i).HR170; daily_PAH_include(i).Hr170(selectweeks(j))];
            % daily_PAH_include_trun(i).weekly_raw = struct([]); % Not sure how to include this
            daily_PAH_include_trun(i).end_expdist_mu = [daily_PAH_include_trun(i).end_expdist_mu; daily_PAH_include(i).end_expdist_mu(selectweeks(j))];
            daily_PAH_include_trun(i).freq_sustained = [daily_PAH_include_trun(i).freq_sustained; daily_PAH_include(i).freq_sustained(selectweeks(j))];
            daily_PAH_include_trun(i).int_expdist_mu = [daily_PAH_include_trun(i).int_expdist_mu; daily_PAH_include(i).int_expdist_mu(selectweeks(j))];
            daily_PAH_include_trun(i).weekly_mean0SCpct = [daily_PAH_include_trun(i).weekly_mean0SCpct; daily_PAH_include(i).weekly_mean0SCpct(selectweeks(j))];
            daily_PAH_include_trun(i).AKcompliance = [daily_PAH_include_trun(i).AKcompliance; daily_PAH_include(i).AKcompliance(selectweeks(j))];
            daily_PAH_include_trun(i).maxofftime = [daily_PAH_include_trun(i).maxofftime; daily_PAH_include(i).maxofftime(selectweeks(j))];
            daily_PAH_include_trun(i).ratio6MWT = [daily_PAH_include_trun(i).ratio6MWT; daily_PAH_include(i).ratio6MWT(selectweeks(j))];
        end
    end

    % Convert truncated structure into cell.
    daily_PAH_cell = struct2cell(daily_PAH_include_trun);

    PAH_hr0 = cell2mat((squeeze(daily_PAH_cell(hr0_ind,1,:)))');
    PAH_hr1 = cell2mat((squeeze(daily_PAH_cell(hr1_ind,1,:)))');
    PAH_sc = cell2mat((squeeze(daily_PAH_cell(sc_ind,1,:)))');
    PAH_map = cell2mat((squeeze(daily_PAH_cell(hr_sc_ind,1,:)))');
    PAH_area = cell2mat(squeeze(daily_PAH_cell(area_ind,1,:)));
    PAH_act_pct = cell2mat(squeeze(daily_PAH_cell(act_pct_ind,1,:)));
    PAH_hour_step = cell2mat((squeeze(daily_PAH_cell(hour_step_ind,1,:)))');
    PAH_max_sc_hr = cell2mat((squeeze(daily_PAH_cell(max_sc_hr_ind,1,:)))'); % excluded
    PAH_peak_pf = cell2mat((squeeze(daily_PAH_cell(peak_pf_ind,1,:)))');
    PAH_MWT = cell2mat((squeeze(daily_PAH_cell(MWT_ind,1,:)))'); % excluded
    PAH_Pred_SC_HR = cell2mat((squeeze(daily_PAH_cell(Pred_SC_HR_ind,1,:)))'); % excluded
    PAH_wear_pct = cell2mat((squeeze(daily_PAH_cell(wear_pct_ind,1,:)))');
    PAH_sustained_walk = cell2mat((squeeze(daily_PAH_cell(sustained_walk_ind,1,:)))');
    PAH_sedentary = cell2mat((squeeze(daily_PAH_cell(sedentary_ind,1,:))));
    PAH_demographics = cell2mat((squeeze(daily_PAH_cell(demographics_ind,1,:)))');
    PAH_demographics = repelem(PAH_demographics,thresh,1);
    PAH_usage_ratio = cell2mat((squeeze(daily_PAH_cell(usage_ratio_ind,1,:)))');

    % Statistical analysis
    M_rf100 = zeros(3,5);
    
    %%% Pared-down version of PCA: Trial 1d (5 variables)
    all_sub_stat = [PAH_hr0(:,1),PAH_hr1(:,2),PAH_sc(:,1:2),PAH_sedentary];

    [max_stat,max_ind] = max(all_sub_stat,[],1);
    [min_stat,min_ind] = min(all_sub_stat,[],1);
    area.min_var = min_stat;
    area.max_var = max_stat;
    
    %%% standardizing variables (mean of 0, SD of 1)
    input_pca2 = normalize(all_sub_stat);
    
    % PCA on "standardized" data
    [coeff_all2,score_all2,latent_all2,~,explained2,~] = pca(input_pca2);
    R2 = corrcoef(input_pca2);
    
    variance_captured(k) = sum(explained2(1:2));
    
end


%% Assessing robustness of PCA

varcap_mean = mean(variance_captured);
varcap_std = std(variance_captured);


%% Find magnitudes of loading vectors

mag = [];
for i = 1:length(coeff_all2)
    mag(i) = sqrt(coeff_all2(i,1)^2 + coeff_all2(i,2)^2);
end
mag = mag';
[maxval, maxind] = maxk(mag,8); % Can change the number of largest elements you want

all_sub_stat = all_sub_stat(:,maxind);
%%% (2) standardizing variables (mean of 0, SD of 1)
input_pca3 = normalize(all_sub_stat);

% PCA on "standardized" data
[coeff_all3,score_all3,latent_all3,~,explained3,~] = pca(input_pca3);
R2 = corrcoef(input_pca3);

%% scree plot of PCA
fh1 = figure('units','normalized','position',[0 0 1 1]);

plot(explained2,'r-','Marker','o','MarkerSize',15,'LineWidth',1.5);

grid on

xlabel('PC ID');
ylabel('Percentage of total variance');

set(gca,'FontSize',15,'FontWeight','bold');

% filename1=[PathName1 '\PCA\ScreenPlot'];
% filename1= 'PCAScreenPlot';
% print(fh1,'-djpeg',filename1);
% close(fh1);
%% plot pca loading overall
fh9 = figure('units','normalized','position',[0 0 1 1]);
b = biplot(coeff_all2(:,1:2),'varlabels',var_name_trunc1d);
xlabel('PC 1','FontSize',15,'FontWeight','bold');
ylabel( 'PC 2','FontSize',15,'FontWeight','bold' );

title('Loading Plot of Variables','FontSize',15,'FontWeight','bold');
set(0,'defaultTextFontSize',15);

% filename9=[PathName1 '\PCA\Loading_area'];
% filename9= 'PCALoadingArea';
% print(fh9,'-djpeg',filename9);
% close(fh9);

%% scatter plot with pca 1&2
group_name = [];
for i = 1:length(daily_PAH_include_trun)
    for j = 1:length(daily_PAH_include_trun(i).hr_mu0)
        group_name = [group_name; convertCharsToStrings(daily_PAH_include_trun(i).name)];
    end
end

glabel = categorical(group_name);
fh10 = figure('units','normalized','position',[0 0 1 1]);

g = gscatter(score_all2(:,1),score_all2(:,2),glabel,'rgbkmc','o+s*d',10,'on');

for k = 1:length(daily_PAH_include_trun)
    g(k).LineWidth = 2;
end

xlabel('PCA1','FontSize',15,'FontWeight','bold');
ylabel( 'PCA2','FontSize',15,'FontWeight','bold' );

set(gca,'FontSize',15,'FontWeight','bold');

%% Write .csv file
writematrix(input_pca, 'input_pca.csv')
writematrix(glabel, 'input_pca_label.csv')
writematrix(input_pca2, 'input_pca_stan.csv')


%% From sparse PCA in Python

sparse_PCA_val = readmatrix('trans_input_pca.csv');
sparse_PCA_com = readmatrix('trans_input_pca_components.csv');
sparse_PCA_val(:,1) = []; sparse_PCA_val(1,:) = [];
sparse_PCA_com(:,1) = []; sparse_PCA_com(1,:) = [];

group_name = [];
for i = 1:length(daily_PAH_include_trun)
    for j = 1:length(daily_PAH_include_trun(i).hr_mu0)
        group_name = [group_name; convertCharsToStrings(daily_PAH_include_trun(i).name)];
    end
end

glabel = categorical(group_name);

% figure(1)
hold on
fh10 = figure('units','normalized','position',[0 0 1 1]);

g = gscatter(sparse_PCA_val(:,1),sparse_PCA_val(:,2),glabel,'rgbkmc','o+s*d',10,'on');

for k = 1:length(daily_PAH_include_trun)
    g(k).LineWidth = 2;
end

xlabel('Sparse PC1','FontSize',15,'FontWeight','bold');
ylabel('Sparse PC2','FontSize',15,'FontWeight','bold' );

set(gca,'FontSize',15,'FontWeight','bold');
hold off


% figure(2)
hold on
fh9 = figure('units','normalized','position',[0 0 1 1]);
b = biplot(sparse_PCA_com(:,1:2),'varlabels',var_name);
xlabel('PC1','FontSize',15,'FontWeight','bold');
ylabel('PC2','FontSize',15,'FontWeight','bold' );

title('Loading Plot of Variables','FontSize',15,'FontWeight','bold');
set(0,'defaultTextFontSize',15);
hold off

