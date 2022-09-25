function [hr_sk0,hr_ku0,hr_ks0,hr_sk1,hr_ku1,hr_ks1,sc_sk,sc_ku,sc_ks] = plot_week_dist2(hr0,hr1,step1,subi_name,Week_No,Year_No)

% make new directory
PathName1 = 'C:\Users\zxu11\OneDrive - Johns Hopkins University\JHU\Projects\Fitbit_IRB_28FEB2018\Plots\';
folder_HR1 = [PathName1 'Weekly Distribution\' 'HR\' subi_name '\'];
folder_HR0 = [PathName1 'Weekly Distribution\' 'HR0\' subi_name '\'];
folder_SC = [PathName1 'Weekly Distribution\' 'SC\' subi_name '\'];

% folder_map = [PathName1 'Weekly Distribution\' 'SC vs HR\' subi_name '\'];

mkdir(folder_HR1)
mkdir(folder_HR0)
mkdir(folder_SC)


%% gaussian distribution function of heart rate at SC=0
pd_hr0 = fitdist(hr0,'Normal');
hr_sk0 = skewness(hr0);
hr_ku0 = kurtosis(hr0);

% measure the goodness of the fit
[~,~,hr_ks0,~] = kstest(hr0,'CDF',pd_hr0);

% plot distribution of heart rate at SC=0
x0 = sort(unique(hr0));

% fh3 = figure('units','normalized','position',[0 0 1 1]);
% histfit(hr0,length(x0),'Normal');
% axis tight
% 
% TextLocation(['Skewness =' num2str(hr_sk0), ', Kurtosis =' num2str(hr_ku0),...
%     ', KS-statistics =' num2str(hr_ks0)],'Location','best');
% 
% xlabel('heart rate','FontSize',15,'FontWeight','bold');
% title( ['Distribution of Heart Rate ' subi_name])
% set(gca,'FontSize',15,'FontWeight','bold');
% filename1=[folder_HR0  subi_name ' Year' num2str(Year_No) ' Week' num2str(Week_No)];
% print(fh3,'-djpeg',filename1);
% close(fh3);

%% gaussian distribution function of heart rate at SC>0
pd_hr1 = fitdist(hr1,'Normal');
hr_sk1 = skewness(hr1);
hr_ku1 = kurtosis(hr1);

% measure the goodness of the fit
[~,~,hr_ks1,~] = kstest(hr1,'CDF',pd_hr1);

% plot distribution of heart rate at SC=1
x1 = sort(unique(hr1));

% fh3 = figure('units','normalized','position',[0 0 1 1]);
% histfit(hr1,length(x1),'Normal');
% axis tight
% 
% TextLocation(['Skewness =' num2str(hr_sk1), ', Kurtosis =' num2str(hr_ku1),...
%     ', KS-statistics =' num2str(hr_ks1)],'Location','best');
% 
% xlabel('heart rate','FontSize',15,'FontWeight','bold');
% title( ['Distribution of Heart Rate ' subi_name])
% set(gca,'FontSize',15,'FontWeight','bold');
% filename1=[folder_HR1  subi_name ' Year' num2str(Year_No) ' Week' num2str(Week_No)];
% print(fh3,'-djpeg',filename1);
% close(fh3);

%% Apply log-normal curve fit to step count distribution without zero step count

pd_step = fitdist(step1,'lognormal');
sc_sk = skewness(step1);

sc_ku = kurtosis(step1);

% measure the goodness of the fit
[~,~,sc_ks,~] = kstest(step1,'CDF',pd_step);

% histogram of step count
x2 = sort(unique(step1));

% fh2 = figure('units','normalized','position',[0 0 1 1]);
% histfit(step1,length(x2),'lognormal');
% axis tight
% TextLocation(['Skewness =' num2str(sc_sk), ', Kurtosis =' num2str(sc_ku),...
%     ', KS-statistics =' num2str(sc_ks)],'Location','best');
% 
% xlabel('step count','FontSize',15,'FontWeight','bold');
% title( ['Distribution of Step Count ' subi_name])
% set(gca,'FontSize',15,'FontWeight','bold');
% 
% filename1 = [folder_SC  subi_name ' Year' num2str(Year_No) ' Week' num2str(Week_No)];
% print(fh2,'-djpeg',filename1);
% close(fh2);

