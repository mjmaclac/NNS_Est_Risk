function[c,ceq]=myconlin(x,use_regions) 
    time_trend2 = repmat((-1:157)',1,length(use_regions));
    c = diff(exp(x(8)+x(end).*time_trend2))-15;
    ceq  = [];
end