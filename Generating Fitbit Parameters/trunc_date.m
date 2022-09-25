function [hr,step,act,t] = trunc_date(n,subi_hr,subi_step,subi_act,subi_time)

% based on each specific date, truncate the datetime and other variables
if n == 4
    start_date = datetime(2018,7,11);
    end_date = datetime(2018,7,11) + days(90);
elseif n == 5
    start_date = datetime(2018,10,15);
    end_date = datetime(2018,12,1);
elseif n == 6
    start_date = datetime(2019,7,9);
    end_date = datetime(2019,7,9) + days(90);
elseif n == 9
    start_date = datetime(2020,3,20);
    end_date = datetime(2021,9,1);
elseif n == 10
    start_date = datetime(2020,3,11);
    end_date = datetime(2021,9,1);
elseif n == 11
    start_date = datetime(2020,3,11);
    end_date = datetime(2021,9,1);
end

start_num = datenum(start_date);
end_num = datenum(end_date);

ind_start = find(subi_time > start_num); 
ind_end = find(subi_time < end_num);

ind = intersect(ind_start,ind_end);

hr = subi_hr(ind);
step = subi_step(ind);
act = subi_act(ind);
t = subi_time(ind);