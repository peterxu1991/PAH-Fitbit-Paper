function [mx_off,num] = plot_fract(val,off)
arr =  val(:);
idx0 =  find(arr == off);
idx0_diff = diff(idx0);
idx1_bigger = find(idx0_diff > 1);

idx1 =find(idx0_diff ==1);

if length(idx1)== 0 &&  length(idx0) > 0
    mx_off = 1;
    num=1;


elseif length(idx1_bigger) == 0
    mx_off = length(idx0);
    num=length(idx0);

 elseif length(idx0) == 0
    mx_off = 0;
    num=0;
else

                                L = zeros(1,length(idx1_bigger)+1);
 
                                for j = 1:length(idx1_bigger)+1
                                    
                                    if j==1
                                        
                                        L(1,j) = idx1_bigger(j);


                                    elseif j <  length(idx1_bigger)+1

                                        L(1,j) = idx1_bigger(j) - idx1_bigger(j-1);
                                
                                    else
                                        
                                        L(1,j) = length(idx0_diff)-idx1_bigger(end) +1;
                                    end
                                
                                end
                                
                                 Y =L;
                                 mx_off = max(Y);
                                 num = length(idx1_bigger)+1;
                               
end
end