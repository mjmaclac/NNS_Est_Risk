function[c,ceq]=mycon(x,Tot_val_ln,Cum_Val_ln,use_regions,len) 
    time_trend2 = repmat((-1:157)',1,length(use_regions));
    c = diff(exp(x(8)+x(23).*time_trend2) ...
                  +repmat(x(9:15)',len,1).*Tot_val_ln(:,use_regions,3) ...
                  +repmat(x(16:22)',len,1).*Cum_Val_ln(:,use_regions,3))-15;
    ceq  = [];
end