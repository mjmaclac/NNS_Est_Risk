%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Confidence intervals and signficane stars %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

parms_mod=parms;   % Copy parameter estimates to be used in constructing CIs
crit_val=1.92;     % Critical value for 95% CI
prec_val=0.001;    % Tolerance for CI search

% Initialize matrices to be populated 
CI_val_lower=zeros(length(parms)+3,1);  
CI_val_upper=zeros(length(parms)+3,1);
LR_check=zeros(length(parms),2);
star_vec=zeros(length(parms)+3,1);

% Use likelihood ratio principles to optimize 
for i=1:length(parms_mod)
    CI_val=eps*sign(parms_mod(i));
    fun_nLL=@(x) nLL_est_pfl_CI(x,a0,Cum_Val_ln,Tot_val_ln,Obs_inv, ...
        effort,tt_ord,use_years,tt_inside_exp,use_regions, ...
        len,len2,t,u,i,CI_val);
    [parms_CI,fval_CI]=fmincon(fun_nLL,parms_mod,A,b,Aeq,beq, ...
        lb', ub',nonlcon,opts); 
    if fval_CI-fval<crit_val && i<16                                                    
            CI_val_lower(i)=0;
    elseif fval_CI-fval<crit_val
        CI_max=0;
        CI_min=-10*parms_mod(i);
        DCI=100;
        while (abs(fval_CI-fval-crit_val)>prec_val && abs(DCI)>abs(parms_mod(i)*prec_val))
            fun_nLL=@(x) nLL_est_pfl_CI(x,a0,Cum_Val_ln,Tot_val_ln, ...
                Obs_inv,effort,tt_ord,use_years,tt_inside_exp, ...
                use_regions,len,len2,t,u,i,CI_val);
            [parms_CI,fval_CI]=fmincon(fun_nLL,parms_mod,A,b,Aeq,beq, ...
                lb', ub',nonlcon,opts); 
            CI_val_lower(i)=CI_val;
            if fval_CI-fval>crit_val
                CI_min=CI_val;
                CI_old=CI_val;
                CI_val=(CI_val+CI_max)/2;
            else
                CI_max=CI_val;
                CI_old=CI_val;
                CI_val=(CI_val+CI_min)/2; 
            end
            DCI=CI_val-CI_old;
            if DCI==0
                CI_max=CI_val;
                CI_min=0;  
                CI_val=(CI_val+CI_min)/2; 
            end
        end
    else
        CI_max=parms_mod(i);
        CI_min=0;
        DCI=100;
        while (abs(fval_CI-fval-crit_val)>prec_val && abs(DCI)>abs(parms_mod(i)*prec_val))
            fun_nLL=@(x) nLL_est_pfl_CI(x,a0,Cum_Val_ln,Tot_val_ln,Obs_inv,effort,tt_ord,use_years,tt_inside_exp,use_regions,len,len2,t,u,i,CI_val);
            [parms_CI,fval_CI]=fmincon(fun_nLL,parms_mod,A,b,Aeq,beq, lb', ub',nonlcon,opts); 
            CI_val_lower(i)=CI_val;
            if fval_CI-fval>crit_val
                CI_min=CI_val;
                CI_old=CI_val;
                CI_val=(CI_val+CI_max)/2;
            else
                CI_max=CI_val;
                CI_old=CI_val;
                CI_val=(CI_val+CI_min)/2; 
            end
            DCI=CI_val-CI_old;
            if DCI==0
                CI_max=CI_val;
                CI_min=0;  
                CI_val=(CI_val+CI_min)/2; 
            end
        end
    end
    LR_check(i,1)=fval_CI-fval;
    fval_CI=0;
    CI_val=2*parms_mod(i);
    CI_max=1e2*parms_mod(i);
    CI_min=parms_mod(i);
    DCI=100;
    cnt_var=0;
    while (abs(fval_CI-fval-crit_val)>prec_val && abs(DCI)>abs(parms_mod(i)*prec_val))
        if (cnt_var==0 && fval_CI-fval<crit_val)
        	fun_nLL=@(x) nLL_est_pfl_CI(x,a0,Cum_Val_ln,Tot_val_ln,Obs_inv,effort,tt_ord,use_years,tt_inside_exp,use_regions,len,len2,t,u,i,CI_val);
            [parms_CI,fval_CI]=fmincon(fun_nLL,parms_mod,A,b,Aeq,beq, lb', ub',nonlcon,opts); 
            CI_val_upper(i)=CI_val;
            CI_val=2*CI_val;
        else
            cnt_var=1;
            fun_nLL=@(x) nLL_est_pfl_CI(x,a0,Cum_Val_ln,Tot_val_ln,Obs_inv,effort,tt_ord,use_years,tt_inside_exp,use_regions,len,len2,t,u,i,CI_val);
            [parms_CI,fval_CI]=fmincon(fun_nLL,parms_mod,A,b,Aeq,beq, lb', ub',nonlcon,opts); 
            CI_val_upper(i)=CI_val;
            if fval_CI-fval>crit_val
                CI_max=CI_val;
                CI_old=CI_val;
                CI_val=(CI_val+CI_min)/2;
            else
                CI_min=CI_val;
                CI_old=CI_val;
                CI_val=(CI_val+CI_max)/2;
            end
            DCI=CI_val-CI_old;
            if DCI==0
                CI_min=CI_val;
                CI_max=1e2*parms_mod(i);  
                CI_val=(CI_val+CI_min)/2; 
            end
        end
    end
    LR_check(i,2)=fval_CI-fval;
end

a_0_fval=squeeze(min(min(fval_array_temp,[],3),[],2));
a_1_fval=squeeze(min(min(fval_array_temp,[],3),[],1));
a_2_fval=squeeze(min(min(fval_array_temp,[],2),[],1))';

a0_CI_min=find(diff(a_0_fval>fval+crit_val)==min(diff(a_0_fval>fval+crit_val)));
a0_CI_max=find(diff(a_0_fval>fval+crit_val)==max(diff(a_0_fval>fval+crit_val)));

CI_val_lower(length(parms)+1)=(a_0_fval(a0_CI_min)*a_0(a0_CI_min)+a_0_fval(a0_CI_min+1)*a_0(a0_CI_min+1))/(a_0_fval(a0_CI_min)+a_0_fval(a0_CI_min+1));
CI_val_upper(length(parms)+1)=(a_0_fval(a0_CI_max)*a_0(a0_CI_max)+a_0_fval(a0_CI_max+1)*a_0(a0_CI_max+1))/(a_0_fval(a0_CI_max)+a_0_fval(a0_CI_max+1));

a1_CI_min=find(diff(a_1_fval>fval+crit_val)==min(diff(a_1_fval>fval+crit_val)));
a1_CI_max=find(diff(a_1_fval>fval+crit_val)==max(diff(a_1_fval>fval+crit_val)));

CI_val_lower(length(parms)+2)=(a_1_fval(a1_CI_min).*a_1(a1_CI_min)+a_1_fval(a1_CI_min+1).*a_1(a1_CI_min+1))./(a_1_fval(a1_CI_min)+a_1_fval(a1_CI_min+1));
CI_val_upper(length(parms)+2)=(a_1_fval(a1_CI_max).*a_1(a1_CI_max)+a_1_fval(a1_CI_max+1).*a_1(a1_CI_max+1))./(a_1_fval(a1_CI_max)+a_1_fval(a1_CI_max+1));

a2_CI_min=min(find(diff(a_2_fval>fval+crit_val)==min(diff(a_2_fval>fval+crit_val))));
a2_CI_max=max(find(diff(a_2_fval>fval+crit_val)==max(diff(a_2_fval>fval+crit_val))));

CI_val_lower(length(parms)+3)=(a_2_fval(a2_CI_min).*a_2(a2_CI_min)+a_2_fval(a2_CI_min+1)*a_2(a2_CI_min+1))/(a_2_fval(a2_CI_min)+a_2_fval(a2_CI_min+1));
CI_val_upper(length(parms)+3)=(a_2_fval(a2_CI_max).*a_2(a2_CI_max)+a_2_fval(a2_CI_max+1)*a_2(a2_CI_max+1))/(a_2_fval(a2_CI_max)+a_2_fval(a2_CI_max+1));

[CI_val_lower, [parms_mod; a0(1) ; a0(2); a0(3)], CI_val_upper]

% Develop vector of stars that relate to the parameters
CI_val_vec = zeros(length(parms),1);
for i=1:length(parms)
    CI_val=eps*sign(parms_mod(i));
    fun_nLL=@(x) nLL_est_pfl_CI(x,a0,Cum_Val_ln,Tot_val_ln,Obs_inv, ...
        effort,tt_ord,use_years,tt_inside_exp,use_regions, ...
        len,len2,t,u,i,CI_val);
    [parms_CI,fval_CI]=fmincon(fun_nLL,parms_mod,A,b,Aeq,beq, ...
        lb', ub',nonlcon,opts); 
    if fval_CI-fval>1.35 && fval_CI-fval<=1.92
        star_vec(i)=1;
    elseif fval_CI-fval>1.92 && fval_CI-fval<=3.32
        star_vec(i)=2;
    elseif fval_CI-fval>3.32
        star_vec(i)=3;
    end
    CI_val_vec(i)=fval_CI;
end

for i=1:3
    a0_CI_temp0=eps*sign(a0(i));
    a0_CI_temp=a0; a0_CI_temp(i)=a0_CI_temp0; 
    try
    fun_nLL=@(x) nLL_est_pfl_CI(parms_mod,a0_CI_temp,Cum_Val_ln, ...
        Tot_val_ln,Obs_inv,effort,tt_ord,use_years,tt_inside_exp, ...
        use_regions,len,len2,t,u,i,CI_val);
    [parms_CI,fval_CI]=fmincon(fun_nLL,parms_mod,A,b,Aeq,beq, ...
        lb',ub',nonlcon,opts); 
    catch
       fval_CI = fval+100;
    end
    if fval_CI-fval>1.35 && fval_CI-fval<=1.92
        star_vec(length(parms)+i)=1;
    elseif fval_CI-fval>1.92 && fval_CI-fval<=3.32
        star_vec(length(parms)+i)=2;
    elseif fval_CI-fval>3.32
        star_vec(length(parms)+i)=3;
    end
end

star_vec