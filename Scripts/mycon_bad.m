function[c,ceq]=mycon_bad(x,Tot_val_ln,Cum_Val_ln,use_regions,len) 
    time_trend2 = repmat((-1:157)',1,length(use_regions));
    c = diff(exp(x(8)+x(15).*time_trend2 ...
                  +repmat([0 x(9) 0 0 x(10) x(11) 0],len,1).*Tot_val_ln(:,use_regions,3) ...
                  +repmat([0 x(12) 0 0 x(13) x(14) 0],len,1).*Cum_Val_ln(:,use_regions,3)))-15;
    ceq  = [];
end