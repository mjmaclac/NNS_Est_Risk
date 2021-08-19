function[nLL]=nLL_est_pfl_CI(x,a0,Cum_Val_ln,Tot_val_ln,Obs_inv,effort,tt_ord,use_years,tt_inside_exp,use_regions,len,len2,t,u,parm_hold,parm_value)
    [d_est,sub_comp] = f_d_est_CI(x,a0,Cum_Val_ln,Tot_val_ln,effort,tt_ord,tt_inside_exp,use_regions,len,len2,t,u,parm_hold,parm_value);
    LL_annual = Obs_inv(1:length(d_est),use_regions).*log(d_est)-d_est;     % Vector of annual log-likelihood values, excluding the constant term.
    nLL = -sum(LL_annual(:));                                               % Pseudo nLL. 
end