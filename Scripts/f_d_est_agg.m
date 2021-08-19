function[d_est, sub_comp]=f_d_est_agg(x,a0,Cum_Val_ln,Tot_val_ln,effort,len,len2,t,u)

    lambda_bars = x(1);
    cons        = x(2);
    Beta        = x(3);
    gamma       = x(4); %(29:35);%(22:28); zeros(7,1); 
    om          = x(end);

    pr_intro_u_disc_t_cond = tril(1./(1+exp(a0(1)-a0(2)*(t-u).^2+a0(3).*repmat(effort(end-len2+1:end),1,len2)))); %CONDITIONAL on not previously disc, Pr intro in year u (column), disc in year t (row), 
    pr_intro_u_NOTdisc_b4t = tril(cumprod(1 - [zeros(1,len2); pr_intro_u_disc_t_cond(1:end-1,:)] ,1));   %Pr intro u, NOT discovered before year t.
    put = pr_intro_u_disc_t_cond .* pr_intro_u_NOTdisc_b4t;  %Pr intro u, disc in t (UNCONDITIONAL)
    time_trend=zeros(159,1); time_trend2 = (-1:157)'; 

    Invas_mean_est=  exp(cons + Beta.*Tot_val_ln(:,3)+gamma.*Cum_Val_ln(:,3)+om.*time_trend2);

    early_intros=repmat(lambda_bars,len2-len+1,1); %_temp
    
    Invas_long     = vertcat(early_intros,Invas_mean_est(2:end,:));
    d_est          = put*Invas_long;               %expected discovery rate for each year
    sub_comp       = {Invas_long, time_trend, pr_intro_u_disc_t_cond, put};  %subcomponents
end