function [act_abs,act_pct, fa, hr_mu0,hr_sig0,hr_sk0,hr_ku0,hr_ks0,hr_mu1,hr_sig1,...
    hr_sk1,hr_ku1,hr_ks1,sc_mu,sc_sig,sc_sk,sc_ku,sc_ks,slope,intercept,area1,...
    max_hr_week,hr_max_sc,max_step_week,sc_max_hr,avg_top_sc,avg_top_hr,...
    avg_6MWT_sc,avg_6MWT_hr,Pred_SC,Pred_HR,weekly_raw]= Fitbit_DataProcess4(subi_hr,subi_step,subi_act,t,subi_name)

%%% Data process

id_week = week(t,'weekofyear');
% name_week = week(t,'shortname');

week_no = unique(week(t,'weekofyear'));

[yr,mo,da] = ymd(t);

no_year = unique(yr);
no_week = [];


for i = 1 : length(no_year)
    
    ind_year = find(yr == no_year(i));
    temp_year = t(ind_year);
    temp_week = unique(week(temp_year,'weekofyear'));
    no_week = [no_week;temp_week];
    year(i).id_week = week(temp_year,'weekofyear');
    
end

%% weekly data analysis

% pre-allocate the matrix
slope = zeros(length(no_week),1);
intercept = zeros(length(no_week),1);

area = zeros(length(no_week),10);
area1 = zeros(length(no_week),10);
area2 = zeros(length(no_week),10);

fa = zeros(length(no_week),1);%fraction of active minutes(SC>0)

fs1 = zeros(length(no_week),1);
fs2 = zeros(length(no_week),1);
fs3 = zeros(length(no_week),1);
fs4 = zeros(length(no_week),1);

act_abs = zeros(length(no_week),4);
act_pct = zeros(length(no_week),4);

total_step = zeros(length(no_week),1);
max_hr_week = zeros(length(no_week),1);
hr_max_sc = zeros(length(no_week),1);
max_step_week = zeros(length(no_week),1);
sc_max_hr = zeros(length(no_week),1);

avg_top_sc = zeros(length(no_week),1);
avg_top_hr = zeros(length(no_week),1);

avg_6MWT_sc = zeros(length(no_week),1);
avg_6MWT_hr = zeros(length(no_week),1);

hr_mu0 = zeros(length(no_week),1);
hr_sig0 = zeros(length(no_week),1);
hr_ks0 = zeros(length(no_week),1);
hr_sk0 = zeros(length(no_week),1);
hr_ku0 = zeros(length(no_week),1);

hr_mu1 = zeros(length(no_week),1);
hr_sig1 = zeros(length(no_week),1);
hr_ks1 = zeros(length(no_week),1);
hr_sk1 = zeros(length(no_week),1);
hr_ku1 = zeros(length(no_week),1);

sc_mu = zeros(length(no_week),1);
sc_sig = zeros(length(no_week),1);
sc_ks = zeros(length(no_week),1);
sc_sk = zeros(length(no_week),1);
sc_ku = zeros(length(no_week),1);

Pred_SC = zeros(length(no_week),1);
Pred_HR = zeros(length(no_week),1);

for k = 1 : length(no_week)
    
    if length(no_year) == 2
        
        year1 = year(1).id_week; year_week1 = unique(year1);
        year2 = year(2).id_week; year_week2 = unique(year2);
        
        if k <= length(year_week1)
            
            ind_week = find(year1 == no_week(k));
            
        elseif k <= length(year_week2) + length(year_week1) && k > length(year_week1)
            
            ind_week = find(year2 == no_week(k)) + length(year1);
            
        end
        
    elseif length(no_year) == 3
        
        year1 = year(1).id_week; year_week1 = unique(year1);
        year2 = year(2).id_week; year_week2 = unique(year2);
        year3 = year(3).id_week; year_week3 = unique(year3);
        
        if k <= length(year_week1)
            
            ind_week = find(year1 == no_week(k));
            
        elseif k <= length(year_week2) + length(year_week1) && k > length(year_week1)
            
            ind_week = find(year2 == no_week(k)) + length(year1);
            
        elseif k <= length(year_week3)+length(year_week2) + length(year_week1) && k > length(year_week2) + length(year_week1)
            
            ind_week = find(year3 == no_week(k)) + length(year1)+ length(year2);
            
        end
        
    elseif length(no_year) == 1
        
        year1 = year(1).id_week;
        ind_week = find(year1 == no_week(k));
        
    end
    
    day_ind = weekday(t(ind_week));
    num_day = unique(day_ind);
    
    if (k == 1 && range(num_day) < 7) || (k == length(no_week) && range(num_day) < 7) || length(day_ind) < 2 * 24 * 60% every day the subject should wear at least 24 hours in one week
        
        continue;
        
    else
        
        % pre-process the missing data
        %         [filled_hr, filled_step, filled_act] = FillMissingFitbitHourlyData(subi_hr(ind),subi_step(ind),subi_act(ind),t(ind));
        
        filled_hr = subi_hr(ind_week);
        filled_step = subi_step(ind_week);
        filled_act = subi_act(ind_week);
        
        %         [day_pct(k), night_pct(k),missing_pct(k)] = det_wear_pct1(t);
        
        total_step(k) = sum(filled_step);
        
        act1 = length(find(filled_act == 0));
        act2 = length(find(filled_act == 1));
        act3 = length(find(filled_act == 2));
        act4 = length(find(filled_act == 3));
        
        act_abs(k,:) = [act1,act2,act3,act4];
        act_pct(k,:) = [act1,act2,act3,act4]/(act1+act2+act3+act4);
        
        %         % max and min heart rate during the week
        %         [max_hr(k),max_hr_ind] = max(filled_hr);
        %         [min_hr(k),min_hr_ind] = min(filled_hr);
        
        %% Analyze data at bin size = 1
        
        step_week = filled_step;
        hr_week = filled_hr;
        
        
        % find the active minutes when step count > 0
        la = length(find(step_week>0));
        fa(k) = la/length(hr_week);
        
        %         fw = length(hr_week)/10080;% fraction of wear time over the week
        %         active_ratio(k) = fa(k)/fw;
        %% peak performance (top 10 step count) & 6mins window around max step
        [max_step_week(k), max_sc_ind] = max(step_week);
        hr_max_sc(k) = hr_week(max_sc_ind);
        
        [max_hr_week(k),max_hr_ind] = max(hr_week);
        sc_max_hr(k) = step_week(max_hr_ind);
        
        [sort_sc,sort_ind] = sort(step_week,'descend');
        sort_hr = hr_week(sort_ind);
        
        % select top 10 step count and get their mean SC&HR
        top10_sc = sort_sc(1:10);avg_top_sc(k) = mean(top10_sc);
        top10_hr = sort_hr(1:10);avg_top_hr(k) = mean(top10_hr);
        
        % detect the continuous 6 mins around max step
        time_week = t(ind_week);
        time_week(end+1) = time_week(end) + hours(10000);% set end point to pick up last consecutive time
        
        temp_sc = step_week;% use subs to pick the SC in consecutive time frame
        temp_hr = hr_week;
        temp_sc(end+1) = 10000;
        temp_hr(end+1) = 10000;
        
        [~,delta_t,~] = hms(diff(time_week));
        end_I1 = find(delta_t ~= 1);% find the breakpoint without consecutive numbers using difference ~= 1       
        start_point = 1;
        MWT_sc = [];
        MWT_hr = [];
        
        for j = 1 : length(end_I1)
            
            end_ind = end_I1(j);
            int_sc = temp_sc(start_point:end_ind); %SC interval with consecutive time
            int_hr = temp_hr(start_point:end_ind); %HR interval with consecutive time
            
            start_point = end_ind + 1;
            
            if length(int_sc) < 6 || length(int_hr) < 6 
                continue;
            else
                
             mean_scj = movmean(int_sc,6);%average 6MWT SC at jth consecutive interval time
             MWT_sc = [MWT_sc;mean_scj];
             mean_hrj = movmean(int_hr,6);%average 6MWT hr at jth consecutive interval time
             MWT_hr = [MWT_hr;mean_hrj];
             
            end
          
        end
        
        if ~isempty(MWT_sc) || ~isempty(MWT_hr)
            
            avg_6MWT_sc(k) = max(MWT_sc);
            avg_6MWT_hr(k) = max(MWT_hr);
            
        else
            warning('no consecutive 6 minutes steps');
        end
              
        %% separate the HR at SC=0 and SC>0
        % determine heart rate at zero step count
        hr0 = hr_week(find(step_week == 0));
        
        hr_mu0(k) = mean(hr0);%mean of HR at SC=0
        hr_sig0(k) = std(hr0);%SD of HR at SC=0
        
        % determine heart rate without zero step count
        hr1 = hr_week(find(step_week >0));
        
        hr_mu1(k) = mean(hr1); %mean of HR at SC>0
        hr_sig1(k) = std(hr1);%SD of HR at SC>0
        
        %% find step count > 0
        step1 = step_week(find(step_week >0));
        
        sc_mu(k) = mean(step1);
        sc_sig(k) = std(step1);
        
        %% determine the time at each step count level:0~30,31~94,95~135,>=136
        
        %         hr_s4 = hr_week(find(step_week >135)); sc_s4 = step_week(find(step_week >135));
        
        %         fs1(k) = length(hr_s1);fs2(k) = length(hr_s2);
        %         fs3(k) = length(hr_s3);fs4(k) = length(hr_s4);
        %% Predict the HR at SC = 170spm using SC = 80+/-5, 100+/-5, 120+/-5 bpm; refer to method of PWC170 (Johnson,BG_1979)       
        hr_s1 = hr_week(find(step_week <85 & step_week >=75)); sc_s1 = step_week(find(step_week <85 & step_week >=75));
        hr_s2 = hr_week(find(step_week <95 & step_week >=85)); sc_s2 = step_week(find(step_week <95 & step_week >=85));
        hr_s3 = hr_week(find(step_week <105 & step_week >=95)); sc_s3 = step_week(find(step_week <105 & step_week >=95));
        
        % extract the mean step count at each interval
        avg_sc1 = mean(sc_s1);std_sc1 = std(sc_s1);
        avg_sc2 = mean(sc_s2);std_sc2 = std(sc_s2);
        avg_sc3 = mean(sc_s3);std_sc3 = std(sc_s3);
        
        avg_hr1 = mean(hr_s1);std_hr1 = std(hr_s1);
        avg_hr2 = mean(hr_s2);std_hr2 = std(hr_s2);
        avg_hr3 = mean(hr_s3);std_hr3 = std(hr_s3);
        
        hr_y = [avg_hr1, avg_hr2,avg_hr3];
        sc_x =[avg_sc1, avg_sc2,avg_sc3];
        
        if isempty(find(isnan(hr_y)==1)) && isempty(find(isnan(sc_x) ==1))
            
              p1 = polyfit(sc_x,hr_y,1);
            Pred_HR(k) = polyval(p1,170);%predict HR at SC=170
            
        end
       %% Predict the step count at HR = 170bpm using HR = 80+/-5, 100+/-5, 120+/-5 bpm; refer to method of PWC170 (Johnson,BG_1979)
       hr_h1 = hr_week(find(hr_week <85 & hr_week >=75)); sc_h1 = step_week(find(hr_week <85 & hr_week >=75));
       hr_h2 = hr_week(find(hr_week <95 & hr_week >=85)); sc_h2 = step_week(find(hr_week <95 & hr_week >=85));
       hr_h3 = hr_week(find(hr_week <105 & hr_week >=95)); sc_h3 = step_week(find(hr_week <105 & hr_week >=95));
        
        % extract the mean step count at each interval
        avg_sc4 = mean(sc_h1);std_sc4 = std(sc_h1);
        avg_sc5 = mean(sc_h2);std_sc5 = std(sc_h2);
        avg_sc6 = mean(sc_h3);std_sc6 = std(sc_h3);
        
        avg_hr4 = mean(hr_h1);std_hr4 = std(hr_h1);
        avg_hr5 = mean(hr_h2);std_hr5 = std(hr_h2);
        avg_hr6 = mean(hr_h3);std_hr6 = std(hr_h3);
        
        hr_x = [avg_hr4, avg_hr5,avg_hr6];
        sc_y =[avg_sc4, avg_sc5,avg_sc6];
        
        if isempty(find(isnan(hr_x)==1)) && isempty(find(isnan(sc_y) ==1))
            
            p2 = polyfit(hr_x,sc_y,1);
            Pred_SC(k) = polyval(p2,170);%predict HR at SC=170      
            
        end 
           
            
        %% calculate the statistical variables of step count and heart rate distribution fit curve
        [temp_year_no,~,~] = ymd(t(ind_week));
        Year_No = unique(temp_year_no);
        Week_No = no_week(k);
        
         weekly_raw(k).SC = step1;
         weekly_raw(k).HR0 = hr0;
         weekly_raw(k).HR = hr1;
         weekly_raw(k).YearNo = Year_No;
         weekly_raw(k).WeekNo = Week_No;
         
         weekly_raw(k).all_SC = step_week;
         weekly_raw(k).all_HR = hr_week;
         weekly_raw(k).all_actlv = filled_act;
         weekly_raw(k).all_time = t(ind_week);
        
        [hr_sk0(k),hr_ku0(k),hr_ks0(k),hr_sk1(k),hr_ku1(k),hr_ks1(k),...
        sc_sk(k),sc_ku(k),sc_ks(k)] = plot_week_dist2(hr0,hr1,step1,subi_name,no_week(k),Year_No);
        
        %% moving average of local maxima and minima to calculate the area of the map
        if length(step_week) == length(hr0)
            
            area(k,m) = 0;
            area1(k,m) = 0;
            area2(k,m) = 0;
            
        elseif length(step_week) < length(hr0) + 10 % if subjects move less than 10 minutes in a week, then exclude the data
            
            continue;
            
        else
            
            %% polyfit of step count vs heart rate
            p = polyfit(step1,hr1,1);
            x = 1 : 0.5: max(step1);
            y = polyval(p,x);
            
            if abs(p(1)) > 5
                
                slope(k) = inf;
                
            else
                
                slope(k) = p(1);
                
            end
            
            intercept(k) = p(2);
            
            % bin size 1
            edges1 = 0 : 1 : max(step_week);
            step_bin_ind1 = discretize(step_week,edges1);
            no_bin1 = unique(step_bin_ind1);
            
            min_hr = zeros(length(no_bin1),1);
            max_hr = zeros(length(no_bin1),1);
            pct_low = zeros(length(no_bin1),1);
            pct_high = zeros(length(no_bin1),1);
            
            % find local min and max at step bin size 1
            
            for i = 1 : length(no_bin1)
                
                bin_ind1 = find(step_bin_ind1 == no_bin1(i));
                %                         bin_step = step_M(bin_ind1);
                bin_hr = hr_week(bin_ind1);
                
                if isempty(bin_hr)
                    
                    warning('No step bin size exist');
                    min_hr(i) = NaN;
                    max_hr(i) = NaN;
                    pct_low(i) = NaN;
                    pct_high(i) = NaN;
                    
                else
                    
                    
                    min_hr(i) = min(bin_hr);
                    max_hr(i) = max(bin_hr);
                    
                    % use 95% percentile data area of each bin to calculate the area of the map
                    pct_high(i) = prctile(bin_hr,97.5);
                    pct_low(i) = prctile(bin_hr,2.5);
                    
                end
            end
            
            bin1 = 1:10;
            
            for m = 1 : length(bin1)
                
                mov_min_hr = movmean(min_hr,bin1(m),'omitnan');
                mov_max_hr = movmean(max_hr,bin1(m),'omitnan');
                
                mov_pct_high = movmean(pct_high,bin1(m),'omitnan');
                mov_pct_low = movmean(pct_low,bin1(m),'omitnan');
                
                mov_step = movmean(no_bin1,bin1(m),'omitnan');
                
                xbdy = [mov_step;flip(mov_step)];
                ybdy1 = [mov_min_hr;flip(mov_max_hr)];
                ybdy2 = [mov_pct_low; flip(mov_pct_high)];
                
                area(k,m) = polyarea(xbdy,ybdy1);
                area2(k,m) = polyarea(xbdy,ybdy2);
                
                % polyfit the two lines on the boundary
                
                p_up = polyfit(mov_step,mov_max_hr,1);
                interp_up = p_up(2);
                y_up = polyval(p_up,x);
                
                p_down = polyfit(mov_step,mov_min_hr,1);
                interp_down = p_down(2);
                y_down = polyval(p_down,x);
                
                inter_point = InterX([x;y_up],[x;y_down]);
                
                if isempty(inter_point)
                    
                    warning('There is no intersection point');
                    yup_max = polyval(p_up,max(step_week));
                    ydown_max = polyval(p_down,max(step_week));
                    
                    area1(k,m) = ((interp_up - interp_down)+(yup_max - ydown_max))*max(step_week)/2;
                    
                else
                    
                    area1(k,m) = (interp_up - interp_down)*inter_point(1,:)/2;
                    
                end
            end
        end
%                         % make new directory
%                         PathName1 = 'C:\Users\zxu11\OneDrive - Johns Hopkins University\JHU\Projects\Fitbit_IRB_28FEB2018\Plots\Area\';
%                         new_folder = [PathName1 subi_name '_bin' num2str(bin1(m))];
%                         mkdir(new_folder)
%         
%                         fh1 = figure('units','normalized','position',[0 0 1 1]);
%                         grid minor
%                         scatter(step_M,hr_M,'MarkerEdgeColor','b','LineWidth',2)
%                         hold on
%                         plot(xbdy, ybdy1, '*-r','LineWidth',2)
%                         plot(xbdy, ybdy2, 'd-k','LineWidth',1)
%                         %plot(max_step,hr_M(ind_max_step),'dr','MarkerSize',10,'LineWidth',2)
%         
%                                     plot(x,y,'k','LineWidth',2)
%                         %             plot(mov_step,mov_max_hr,'dg','MarkerSize',10,'LineWidth',2)
%                         plot(x,y_up,'--m','LineWidth',2)
%                         %             plot(mov_step,mov_min_hr,'dr','MarkerSize',10,'LineWidth',2)
%                         plot(x,y_down,'-.g','LineWidth',2)
%                         legend('SC_HR_Map','min-max boundary','outlier95% boundary','upper line','lower line','Interpreter', 'none','Location','Best')
%                         grid on
%                         %             s = sprintf('y = (%.2f) x + (%.2f)',p(1),p(2));
%                         %             TextLocation(s,'Location','southeast')
%                         %             text(x(end),y(1),['y = ' num2str(round(p(1),2)) 'x + ' num2str(round(p(2),2))],'FontSize',20)
%                         xlabel('Step count','FontSize',25,'FontWeight','bold')
%                         ylabel('Heart rate','FontSize',25,'FontWeight','bold')
%                         title(['Step & HR at week' num2str(k-1 + min(no))],'FontSize',25,'FontWeight','bold')
%                         set(gca,'FontSize',12,'FontWeight','bold');
%                         hold off
%                         filename1=[new_folder '\'  'week' num2str(k-1 + min(no))];
%                         print(fh1,'-djpeg',filename1);
%                         close(fh1);
        
    end
end
%% remove the zero elment
[ind0,~] = find(area1(:,10) ==0);
[ind1,~] = find(slope == inf);
comb_ind = [ind0;ind1];

fa(comb_ind) = [];

%     fs1(comb_ind) = [];
%     fs2(comb_ind) = [];
%     fs3(comb_ind) = [];
%     fs4(comb_ind) = [];

slope(comb_ind,:) = [];
intercept(comb_ind,:) = [];

area1(comb_ind,:) = [];
act_abs(comb_ind,:) = [];
act_pct(comb_ind,:) = [];

hr_mu0(comb_ind,:) = [];
hr_sig0(comb_ind,:) = [];
hr_sk0(comb_ind,:) = [];
hr_ku0(comb_ind,:) = [];
hr_ks0(comb_ind,:) = [];

hr_mu1(comb_ind,:) = [];
hr_sig1(comb_ind,:) = [];
hr_sk1(comb_ind,:) = [];
hr_ku1(comb_ind,:) = [];
hr_ks1(comb_ind,:) = [];

sc_mu(comb_ind,:) = [];
sc_sig(comb_ind,:) = [];
sc_sk(comb_ind,:) = [];
sc_ku(comb_ind,:) = [];
sc_ks(comb_ind,:) = [];

max_hr_week(comb_ind,:) = [];
hr_max_sc (comb_ind,:) = [];
max_step_week(comb_ind,:) = [];
sc_max_hr(comb_ind,:) = [];

avg_top_sc(comb_ind,:) = [];
avg_top_hr(comb_ind,:) = [];

avg_6MWT_sc(comb_ind,:) = [];
avg_6MWT_hr(comb_ind,:) = [];

Pred_SC(comb_ind,:) = [];
Pred_HR(comb_ind,:) = [];
%% plot the sample entropy of step and heart rate
% PathName1 = 'C:\Users\zxu11\OneDrive - Johns Hopkins University\JHU\Projects\Fitbit_IRB_28FEB2018\Plots\SampleEntropy\';
% fh1 = figure('units','normalized','position',[0 0 1 1]);
% grid minor
% plot(step_en,'-dr','MarkerSize',10,'LineWidth',2)
% hold on
% plot(hr_en,'-og','MarkerSize',10,'LineWidth',2)
% title('Sample Entropy of Weekly Raw Signal')
% legend('Step Count','Heart Rate','Interpreter', 'none','Location','Best')
% xlabel('Week Number','FontSize',25,'FontWeight','bold')
% ylabel('Sample Entropy','FontSize',25,'FontWeight','bold')
% set(gca,'FontSize',12,'FontWeight','bold');
% hold off
% filename1=[PathName1 subi_name ];
% print(fh1,'-djpeg',filename1);
% close(fh1);