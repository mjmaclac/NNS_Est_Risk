function[nLL]=nLL_est_pfl_agg(x,a0,Cum_Val_ln,Tot_val_ln,Obs_inv,effort,len,len2,t,u,use_regions)
    [d_est,sub_comp] = f_d_est_agg(x,a0,Cum_Val_ln,Tot_val_ln,effort,len,len2,t,u);
    LL_annual = sum(Obs_inv((end-212):end,use_regions),2).*log(d_est)-d_est;     % Vector of annual log-likelihood values, excluding the constant term.
    nLL = -sum(LL_annual(:));                                                    % Pseudo nLL. 
end