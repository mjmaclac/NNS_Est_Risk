function[d_est, sub_comp]=f_d_est_CI(x,a0,Cum_Val_ln,Tot_val_ln,effort,tt_ord,tt_inside_exp,use_regions,len,len2,t,u,parm_hold,parm_value)

    x(parm_hold)=parm_value;

    lambda_bars = x(1:7);
    cons        = x(8);
    Beta        = x(9:15);
    gamma       = x(16:22);%(29:35);%(22:28);
    om          = x(end);
    
%     pr_intro_u_disc_t_cond = tril(1./(1+exp(a0(1)-pi*(a0(2)*(t-u)).^2+a0(3).*repmat(effort(end-len2+1:end),1,len2)))); %CONDITIONAL on not previously disc, Pr intro in year u (column), disc in year t (row), 
        pr_intro_u_disc_t_cond = tril(1./(1+exp(a0(1)-a0(2)*(t-u).^2+a0(3).*repmat(effort(end-len2+1:end),1,len2)))); %CONDITIONAL on not previously disc, Pr intro in year u (column), disc in year t (row), 
    pr_intro_u_NOTdisc_b4t = tril(cumprod(1 - [zeros(1,len2); pr_intro_u_disc_t_cond(1:end-1,:)] ,1));   %Pr intro u, NOT discovered before year t.
    put = pr_intro_u_disc_t_cond .* pr_intro_u_NOTdisc_b4t;  %Pr intro u, disc in t (UNCONDITIONAL)
    
    time_trend=zeros(159,1);
    time_trend2 = repmat((1:159)',1,length(use_regions)); 
    Invas_mean_est = exp(cons ...
                            + repmat(Beta',len,1).*Tot_val_ln(:,use_regions,3) ...
                            + repmat(gamma',len,1).*Cum_Val_ln(:,use_regions,3) ...
                            + om.*time_trend2);
    
    early_intros=repmat(lambda_bars',len2-len+1,1); %_temp
        
    Invas_long     = vertcat(early_intros,Invas_mean_est(2:end,:));
    d_est          = put*Invas_long;               %expected discovery rate for each year
    sub_comp       = {Invas_long, time_trend, pr_intro_u_disc_t_cond, put,Invas_long};  %subcomponents
end