clear;clc;close all;
%% read data
PathName = 'C:\Users\zxu11\OneDrive - Johns Hopkins University\JHU\Projects\Fitbit_IRB_28FEB2018\';
FileName = 'PC_Step_HR_11012020.csv';

% Read the csv file with header
txt = readtimetable([PathName FileName]);

%% Detect the key variables we need.

uid = txt.fitbit_uid;
ts = txt.timestamp;%time stamp
hr = txt.heart_rate;%heart rate
sc = txt.step_count;%step count
act = txt.activity_level;%activity level

[~,idx] = sort(datenum(ts),'ascend');
Data_all = [hr,sc,act];
Data_sorted = Data_all(idx,:);

uid_sorted = uid(idx);
time_sorted = ts(idx);
%% Obtain uid list and its corresponding name
% if contains(FileName(1:3),'_')
%    gtype = FileName(1:2);
% else
%    gtype = FileName(1:3);
% end
[num_users,user_txt,~] = xlsread('PC_Users3.csv');
[~,ind6]=find(strcmp(user_txt,'fitbit_uid'));
[~,ind7]=find(strcmp(user_txt,'notes'));
[~,ind8]=find(strcmp(user_txt,'age'));

uid_name_list = user_txt(:,[ind6,ind7]);
user_age = num_users(:,ind8);
%% combine the data with same fitbit uid and then sort them with time stamp

for j = 2 : size(uid_name_list,1)
    
    if isempty(uid_name_list{j,1})
        continue
    end
    x = find(strcmp(uid_sorted,uid_name_list{j,1}));
    
    if isempty(x)
        warning('The user is not on the list');
        continue
    else
                
        RawFBData(j-1).step = Data_sorted(x,2);
        RawFBData(j-1).hr = Data_sorted(x,1);
        RawFBData(j-1).act = Data_sorted(x,3);
        RawFBData(j-1).name = uid_name_list{j,2};
        RawFBData(j-1).age = user_age(j-1);
        
        for i = 1 : length(x)
            RawFBData(j-1).time(i,1) = datenum(time_sorted(x(i)));
        end
    end
    
end

save 'PC_Step_HR_110220.mat' RawFBData;