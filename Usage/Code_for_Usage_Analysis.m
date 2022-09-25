clc;clear all;close all;
period =7;
idx =[15,17,19,20,21,22,23,26,27,28,30];



Final =[];
for i=1:length(idx)
    if idx(i)<10
        data = xlsread('FBPH0'+string(idx(i))+'.xlsx');
    else
        data = xlsread('FBPH'+string(idx(i))+'.xlsx');
    end
    
    
    cap = floor(size(data,2)/period);
    
    
    
    FF =[];
    
    for j =1:cap
            if j==cap
                val = data(2:end,(j-1)*period+1:end);
            else
                val = data(2:end,(j-1)*period+1:j*period);
            end

            
            [a0,kk] = plot_fract(val,0);
            
            b0 = mean(val(:));

            FF = [FF;[a0,b0]];          
    end

   
    B = FF(:,2);
    A = FF(:,1);
    p_A = polyfit( 1:length(A) , A , 1);
    B_pr  =  ewn(B,0.3);
    p_B = polyfit( 1:length(B) , B_pr , 1);
    
    array = [mean(A),std(A),p_A(1),mean(B),std(B),p_B(1),idx(i)];
    Final = [Final;array];
   

end


col_header = {'Avergae Max.Off time','STD Max.Off time','Slope Max.Off time','Average Compliance','STD Compliance','Slope Compliance','Subject'};
Final(:,1:5) = round(Final(:,1:5),3);
A0 = num2cell(Final);
A0 = [col_header;A0];
xlswrite('Usage_Analysis.xlsx',A0)