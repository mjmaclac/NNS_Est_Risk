function[d_est, sub_comp]=f_d_est(x,a0,Cum_Val_ln,Tot_val_ln,effort,tt_ord,tt_inside_exp,use_regions,len,len2,t,u)

    lambda_bars = x(1:7);
    cons        = x(8);
    if length(x)>9
        Beta        = x(9:15);
    end
    if length(x)>9
        gamma       = x(16:22); %(29:35);%(22:28); zeros(7,1); 
    end
%     Beta    = zeros(7,1);
%     gamma   = Beta; %(29:35);%(22:28); zeros(7,1); 
%     Beta([2 5 6])    = x(9:11);
%     gamma([2 5 6])   = x(12:14); %(29:35);%(22:28); zeros(7,1); 
    if length(x)==24
        om2         = x(end-1);
    end
    om          = x(end);

    pr_intro_u_disc_t_cond = tril(1./(1+exp(a0(1)-a0(2)*(t-u).^2+a0(3).*repmat(effort(end-len2+1:end),1,len2)))); % CONDITIONAL on not previously disc, Pr intro in year u (column), disc in year t (row), 
    pr_intro_u_NOTdisc_b4t = tril(cumprod(1 - [zeros(1,len2); pr_intro_u_disc_t_cond(1:end-1,:)] ,1));            % Pr intro u, NOT discovered before year t.
    put = pr_intro_u_disc_t_cond .* pr_intro_u_NOTdisc_b4t;                                                       % Pr intro u, disc in t (UNCONDITIONAL)
    time_trend=zeros(159,1);
    time_trend2 = repmat((-1:157)',1,length(use_regions)); 
    time_trend3 = repmat([zeros(68,1); (1:91)'],1,length(use_regions));

    if length(x)==23
        Invas_mean_est=  exp(cons ...
                            +repmat(Beta',len,1).*Tot_val_ln(:,use_regions,3) ...
                            +repmat(gamma',len,1).*Cum_Val_ln(:,use_regions,3) ...
                            +om.*time_trend2);
    elseif length(x)==24
        Invas_mean_est=  exp(cons ...
                            +repmat(Beta',len,1).*Tot_val_ln(:,use_regions,3) ...
                            +repmat(gamma',len,1).*Cum_Val_ln(:,use_regions,3) ...
                            +om.*time_trend2 + om2.*time_trend3);        
    elseif length(x)==9
        Invas_mean_est = exp(cons +om.*time_trend2);
    else
        Invas_mean_est=  exp(cons ...
                +repmat(Beta',len,1).*Tot_val_ln(:,use_regions,3) ...
                +om.*time_trend2);
    end

    early_intros = repmat(lambda_bars',len2-len+1,1); 
    Invas_long   = vertcat(early_intros,Invas_mean_est(2:end,:));

    d_est          = put*Invas_long;               %expected discovery rate for each year
    sub_comp       = {Invas_long, time_trend, pr_intro_u_disc_t_cond, put};  %subcomponents
end