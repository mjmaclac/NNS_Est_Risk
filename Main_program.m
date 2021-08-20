%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Model of trade, introduction and discovery %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear all; close all;                     % Clear stored information and workspace
dir=fileparts(which('Main_program'));          % Set to location of main program

% Paths to main directory and all sub-directories 
cd(dir); 
export_dir=[dir '\Figures\'];      scripts_dir=[dir '\Scripts\'];  
spread_dir=[dir '\Spreadsheets\']; output_dir=[dir '\Output\']; 
addpath(genpath(dir)) 

%% Bring in data and set globals 
try 
    load([output_dir 'import_inv_data.mat']);   % Load organized data if already available
catch                                           % If not, organize data and save
    Fit_effort                                  % Construct proxy for effort 
    Import_values                               % Trade information
    Agg_setup                                   % Aggregate model setup (including invasive species data)
    save([output_dir 'import_inv_data'])        % Save information
end

%% Regionally aggregated model
try 
    load([output_dir 'Agg_results.mat'])
    pr_intro_u_disc_t_cond_array=zeros(len2,len2,profs_tot);            % ... profile probability of observation
    put_array=zeros(len2,len2,profs_tot);                               % ... probability of discovery
catch
    t=repmat([1:len2]',1,len2); u=repmat([1:len2],len2,1);              % Define time indices
    tt_ord=2;                                                           % Order of spline (not currently applied)
    effort=effort((end-len2+1):end);                                    % Truncate effort variable 
    Tot_val_inc=squeeze(Tot_val_ln(:,1,:));                             % Annual imports 
    Cum_Val_inc=squeeze(Cum_Val_ln(:,1,:));                             % Cumulative imports

    %     kappa c,   beta gamma omega 
    x0 = [7.2   3.5  0.5  -0.1  0]';                                    % Initial guess for parameters
    lb = [0    -1000 0    -1    -100]';                                 % Lower bound
    ub = [100  1000  1    0     100]';                                  % Upper bound
    A=[]; b=[]; nonlcon=[];                                             % (Non-)linear constraints

    profs_num=25;  %25 max,used in paper                                % Resolution of profile--evenly spaced
    alpha_dims=3;                                                       % Number of Delay2Discovery parameters
    profs_tot=profs_num^alpha_dims;                                     % Combinations of alpha parameters
    parms_array=zeros(length(lb),profs_tot);                            % Initialize array of parameter estimates 
    fval_array=zeros(profs_tot,1);                                      % ... pseudo-nLL 
    Invas_mean_est_array=zeros(len2,profs_tot);                         % ... estimated establishments 
    time_trend_array=zeros(len,profs_tot);                              % ... estimated time-trend 
    pr_intro_u_disc_t_cond_array=zeros(len2,len2,profs_tot);            % ... profile probability of observation
    put_array=zeros(len2,len2,profs_tot);                               % ... probability of discovery

    % Vectors of possible values of a0 
    a_0=linspace(1,10,profs_num);           
    a_1=linspace(8.5530e-04,0.002,profs_num); 
    a_2=linspace(-0.02,0,profs_num);     

    opts = optimoptions(@fmincon,'MaxFunEvals',10000, ...
        'MaxIterations',40000,'TolX',1e-30, ...
        'TolFun',1e-30,'Display','off');                                % Options for fmincon

    WaitMessage=parfor_wait(profs_tot,'Waitbar',true);
    
    parfor m=1:profs_tot                                                % Parallel loop over possible combos of alpha
        [i,j,k]=ind2sub(profs_num*ones(1,alpha_dims),m);                % Get indexes out of m
        a0=[a_0(i) a_1(j) a_2(k)];         
        fun_nLL=@(x) nLL_est_pfl_agg(x,a0,Cum_Val_inc,Tot_val_inc, ...  
         Obs_inv,effort,len,len2,t,u,use_regions);                      % Pseudo profile nLL function for given alpha
        [parms,fval,e_flag]=fmincon(fun_nLL,x0,A,b,Aeq,beq,lb',ub', ...
                            nonlcon,opts)
        if isempty(parms)~=1                                            % In case optimization fails...
            [d_est,sub_comp] = f_d_est_agg(parms,a0,Cum_Val_inc, ...
                Tot_val_inc,effort,len,len2,t,u);                       % Recover discovery estimates and other features
            parms_array(:,m)=parms;                                     % Parameters
            fval_array(m)=fval;                                         % Pseudo profile nLL value
            Invas_mean_est_array(:,m)=sub_comp{1};                      % Establishments 
            time_trend_array(:,m)=sub_comp{2};                          % Trend
            pr_intro_u_disc_t_cond_array(:,:,m)=sub_comp{3};            % Probability of observation
            put_array(:,:,m)=sub_comp{4};                               % Probability of discovery
        end
        WaitMessage.Send;
    end

    WaitMessage.Destroy
    
    % Display outputs 
    fval_array(fval_array==0)=1e6;
    fval_array_temp=reshape(fval_array,profs_num*ones(1,alpha_dims));
    [val_temp,ind_temp]=min(fval_array_temp(:));  
    [i,j,k]=ind2sub(profs_num*ones(1,alpha_dims),ind_temp);
    a0=[a_0(i) a_1(j) a_2(k)];
    fun_nLL=@(x) nLL_est_pfl_agg(x,a0,Cum_Val_inc,Tot_val_inc, ...     
         Obs_inv,effort,len,len2,t,u,use_regions);    
    [parms,fval]=fmincon(fun_nLL,x0,A,b,Aeq,beq, lb', ub',nonlcon,opts);   

    [d_est,sub_comp] = f_d_est_agg(parms,a0,Cum_Val_inc, ...
                Tot_val_inc,effort,len,len2,t,u);  
    Invas_mean_est=sub_comp{1}; Invas_long=sub_comp{1};
    time_trend=sub_comp{2};
    pr_intro_u_disc_t_cond=sub_comp{3}; put=sub_comp{4};
    save([output_dir 'Agg_results.mat'])
end

%% Aggregate figures for SI Appendix
Agg_figs_SIA

%% Regionally disaggregated data and setup
Reg_NNs           
    
%% Regionally disaggregated model
try 
    load([output_dir 'Regional_results.mat'])
    pr_intro_u_disc_t_cond_array=zeros(len2,len2,profs_tot);            % Defined below. Too large to save. 
    put_array=zeros(len2,len2,profs_tot);                               % ""
catch
    % Bounds and intial value--before period and slopes
    t=repmat([1:len2]',1,len2); u=repmat([1:len2],len2,1);
    tt_ord=2;
    effort=effort((end-len2+1):end);
    Tot_val_inc=Tot_val_ln; Cum_Val_inc=Cum_Val_ln; 
    
         %kappas        constant betas         gammas           omega
    x0 = [0.1*ones(1,7)  0       0.5*ones(1,7) zeros(1,7)       0]'; 
    lb = [zeros(1,7)    -10000   zeros(1,7)    -10000*ones(1,7) -100]';
    ub = [100*ones(1,7)  10000   10*ones(1,7)  10000*ones(1,7)  100]'; 
    A=[]; b=[];  
    nonlcon = @(parms_late) mycon(parms_late,Tot_val_ln,Cum_Val_ln,...
         use_regions,len);

% Alternative starting values and linear/nonlinear constraints for
% robustness checks described in the manuscript

% Non-trade only
%     x0 = [0.1*ones(1,7) 0 0]';         
%     lb = [zeros(1,7)    -10000     -1000]';
%     ub = [100*ones(1,7)  10000     1000]'; 
%     A=[]; b=[];  
%     nonlcon = @(parms_late) myconlin(parms_late,use_regions);

% Exclude regional parameters when beta and gamma coefficents not
% significant
%     x0 = [0.1*ones(1,7)  0         0.5*ones(1,3) zeros(1,3)      0]';  
%     lb = [zeros(1,7)    -10000     zeros(1,3)    -10000*ones(1,3) -100]';
%     ub = [100*ones(1,7)  10000     10*ones(1,3)  10000*ones(1,3)  100]'; 
%     A=[]; b=[]; 
%     nonlcon = @(parms_late) mycon_bad(parms_late,Tot_val_ln,Cum_Val_ln,use_regions,len);

    profs_num=25;  %~25 max for machines with 16GB RAM                  % Resolution of profile--evenly spaced
    alpha_dims=3;                                                       % Number of discovery parameters
    profs_tot=profs_num^alpha_dims;                                     % Combinations of alpha parameters
    parms_array=zeros(length(lb),profs_tot);                            % Initialize array of parameter estimates 
    fval_array=zeros(profs_tot,1);                                      % ""         pseudo-nLL 
    Invas_mean_est_array=zeros(len2,7,profs_tot);                       % ""         estimated establishments 
    time_trend_array=zeros(len,profs_tot);                              % ""         estimated time-trend 
    pr_intro_u_disc_t_cond_array=zeros(len2,len2,profs_tot);            % ""         profile probability of observation
    put_array=zeros(len2,len2,profs_tot);                               % ""         probability of discovery
    
    a_0=linspace(1,10,profs_num);                                       % Vector of possible values of a0 values 
    if length(x0)>=23                                                   % Vector of possible a1 values---depends on specification
        a_1=linspace(1e-04,0.002,profs_num); 
    else
        a_1=linspace(-2e-04,2e-04,profs_num);
    end
    a_2=linspace(-0.02,0,profs_num);                                    % Vector of possible a2 values

    opts = optimoptions(@fmincon,'MaxFunEvals',10000, ...
        'MaxIterations',40000,'TolX',1e-30, ...
        'TolFun',1e-30,'Display','off');                                % Options for fmincon

    WaitMessage=parfor_wait(profs_tot,'Waitbar',true);                  % Start wait message
    
    parfor m=1:profs_tot                                                % Parallel loop over possible combos of alpha
        [i,j,k]=ind2sub(profs_num*ones(1,alpha_dims),m);                % Get indexes out of m
        a0=[a_0(i) a_1(j) a_2(k)];         
        fun_nLL=@(x) nLL_est_pfl(x,a0,Cum_Val_inc,Tot_val_inc, ...      % divide by 1e12 if using non-logged cumulative trade. 
         Obs_inv,effort,tt_ord,use_years,tt_inside_exp, ...
         use_regions,len,len2,t,u);                                     % Pseudo profile nLL function for given alpha
        [parms,fval,e_flag]=fmincon(fun_nLL,x0,A,b,Aeq,beq,lb',ub', ...
                            nonlcon,opts);
        if isempty(parms)~=1                                            % In case optimization fails...
            [d_est,sub_comp] = f_d_est(parms,a0,Cum_Val_inc, ...
                Tot_val_inc,effort,tt_ord,tt_inside_exp, ...
                use_regions,len,len2,t,u);                              % Recover discovery estimates and other features
            parms_array(:,m)=parms;                                     % Parameters
            fval_array(m)=fval;                                         % Pseudo profile nLL value
            Invas_mean_est_array(:,:,m)=sub_comp{1};                    % Establishments 
            time_trend_array(:,m)=sub_comp{2};                          % Trend
            pr_intro_u_disc_t_cond_array(:,:,m)=sub_comp{3};            % Probability of observation
            put_array(:,:,m)=sub_comp{4};                               % Probability of discovery
        end
        WaitMessage.Send;
    end

    WaitMessage.Destroy

    % Format results    
    fval_array(fval_array==0) = 1e6;
    fval_array_temp           = reshape(fval_array,profs_num*ones(1,alpha_dims));
    [val_temp,ind_temp]       = min(fval_array_temp(:)); 
    [i,j,k]                   = ind2sub(profs_num*ones(1,alpha_dims),ind_temp);
    a0                        = [a_0(i) a_1(j) a_2(k)];
    fun_nLL=@(x) nLL_est_pfl(x,a0,Cum_Val_inc,Tot_val_inc,Obs_inv,effort,tt_ord,use_years,tt_inside_exp,use_regions,len,len2,t,u);
    [parms,fval]              = fmincon(fun_nLL,x0,A,b,Aeq,beq, lb', ub',nonlcon,opts);   

    [d_est,sub_comp] = f_d_est(parms,a0,Cum_Val_inc,Tot_val_inc,effort,tt_ord,tt_inside_exp,use_regions,len,len2,t,u); 
    Invas_mean_est=sub_comp{1};
    time_trend=sub_comp{2};
    pr_intro_u_disc_t_cond=sub_comp{3};
    put=sub_comp{4};
    Invas_long=sub_comp{1};
    save([output_dir 'Regional_results.mat'])                              
end

%% Regional model visulations
Figs_Main

SI_Figures

%% Post-estimation programs   
CI_sign            % Prints point estimates, 95% confidence intervals, and the number of stars 

Delay2Discovery    % Provides the 2.5th median and 97.5th percentiles of the delay to discovery
