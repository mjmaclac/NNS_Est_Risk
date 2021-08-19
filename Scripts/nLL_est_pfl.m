function[nLL]=nLL_est_pfl(x,a0,Cum_Val_ln,Tot_val_ln,Obs_inv,effort,tt_ord,use_years,tt_inside_exp,use_regions,len,len2,t,u)
    [d_est,sub_comp] = f_d_est(x,a0,Cum_Val_ln,Tot_val_ln,effort,tt_ord,tt_inside_exp,use_regions,len,len2,t,u);
    LL_annual = Obs_inv(1:length(d_est),use_regions).*log(d_est)-d_est;     % Vector of annual log-likelihood values, excluding the constant term.
    nLL = -sum(LL_annual(:));                                               % Pseudo nLL. 
    Invas_temp=sub_comp{1};
%     if sum(sum(15<diff(Invas_temp(54:end,:),1)))>0
%         nLL = nLL + 1000;
%     end
end