%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setting up aggregate program  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Time_span     = 1854:2012;              % Vector of years included in sample
tsb_1         = Time_span(1)-1; 
tt_ord        = 2;                      % Time trend order: 2=linear, 3=quad, 4=cubic (not currently used)
tt_inside_exp = 1;                      % Time trend inside exp() or outside
use_years     = (1854:2012) - 1853;     % Years incorporated in likelihood function. Subtract 1853 to convert to count years from actual years
len           = length(use_years);      % Time span when we have import data 1854-2012
len2          = len+54;                 % 1800-2012
A=[]; b=[]; Aeq=[];beq=[]; nonlcon=[];  % Constraints
use_regions=[2:6 8 9];                  % Use regions listed in manuscript

%% Effort 
load 'effort.mat'
effort=effort(end-length(Time_span)-77:end);
effort=effort-mean(effort);

%% Regional Non-native Species 
[num, txt]  = xlsread([spread_dir 'trimmed_reg.xls'], 'July_2021_update');
Obs_inv_reg = zeros(2012-1775,9);

for i=1:length(txt)
    for t=1776:2012
        for r=1:9
            if strcmp(txt(i),reg_list(r)) && num(i,1)==t
                Obs_inv_reg(t-1775,r)=Obs_inv_reg(t-1775,r)+num(i,2);
            end
        end
    end
end

%% Import data
Tot_val_hold           = squeeze(imports_vol_nan(:,:,:));
[all_trade_nums, text] = xlsread([spread_dir 'All_trade'],'Sheet1');
Tot_val                = cat(3,Tot_val_hold, zeros(159,9,1));
for i=1:size(all_trade_nums,1)
    for t=1:size(Tot_val,1)
        for r=2:9
            if strcmp(text(i),reg_list(r)) && all_trade_nums(i,1)==t+1853
                Tot_val(t,r,5)=all_trade_nums(i,3);
            end
        end
    end
end
Tot_val(:,1,5) = sum(Tot_val(:,:,5),2);
[num, txt]     = xlsread([spread_dir 'inflation.xlsx'], 'Sheet4');

for t=(1854:2012)-1853
    for r=2:9
        Tot_val(t,r,5)=Tot_val(t,r,5)*num(t,2);
    end
end
Tot_val(:,:,5)  = max(0,Tot_val(:,:,5)/1000-Tot_val(:,:,3));
Cum_Val_hold    = cumsum(Tot_val,1);
Cum_Val         = cat(1,zeros(1,9,5),Cum_Val_hold(1:end-1,:,:));

Tot_val_ln      = log(max(Tot_val*1e9,1));
Cum_Val_ln_hold = log(cumsum(max(Tot_val*1e9,1),1));
Cum_Val_ln      = cat(1,zeros(1,9,5),Cum_Val_ln_hold(1:end-1,:,:));

Obs_inv         = Obs_inv_reg;
Obs_inv(:,1)    = sum(Obs_inv(:,use_regions),2);
Cum_obs_inv     = cumsum(Obs_inv,1);   