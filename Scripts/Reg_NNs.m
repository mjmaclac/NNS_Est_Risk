%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Import and organize data for estimation of regional model %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

use_regions=[2:6 8 9];
[num, txt]=xlsread([spread_dir 'trimmed_reg.xls'], 'July_2021_update');
Obs_inv_reg=zeros(2012-1800+1,9);

for i=1:length(txt)
    for t=1800:2012
        for r=1:9
            if strcmp(txt(i),reg_list(r)) && num(i,1)==t
                Obs_inv_reg(t-1799,r)=Obs_inv_reg(t-1799,r)+num(i,2);
            end
        end
    end
end

Tot_val_hold=squeeze(imports_vol_nan(:,:,:));
[all_trade_nums, text]=xlsread([spread_dir 'All_trade'],'Sheet1');
Tot_val=cat(3,Tot_val_hold, zeros(159,9,1));
for i=1:size(all_trade_nums,1)
    for t=1:size(Tot_val,1)
        for r=2:9
            if strcmp(text(i),reg_list(r)) && all_trade_nums(i,1)==t+1853
                Tot_val(t,r,5)=all_trade_nums(i,3);
            end
        end
    end
end
Tot_val(:,1,5)=sum(Tot_val(:,:,5),2);
[num, txt]=xlsread([spread_dir 'inflation.xlsx'], 'Sheet4');

for t=(1854:2012)-1853
    for r=2:9
        Tot_val(t,r,5)=Tot_val(t,r,5)*num(t,2);
    end
end
Tot_val                   = Tot_val*1e9;
Tot_val_half              = zeros(size(Tot_val));
Tot_val_half(end,:,:)     = Tot_val(end,:,:);
Tot_val_half(1:end-1,:,:) = (max(Tot_val(1:end-1,:,:),1)+max(Tot_val(2:end,:,:),1))/2;

Cum_Val_hold=cumsum(Tot_val,1);

Tot_val_ln       = log(max(Tot_val,1));
Cum_Val_ln_hold1 = cumsum(Tot_val,1);
Cum_Val_ln_hold2 = max(Cum_Val_ln_hold1,1);
Cum_Val_ln_hold3 = log(Cum_Val_ln_hold2);
Cum_Val_ln       = cat(1,zeros(1,9,5),Cum_Val_ln_hold3(1:end-1,:,:));

Obs_inv         = Obs_inv_reg;
Obs_inv(:,1)    = sum(Obs_inv(:,use_regions),2);
Cum_obs_inv     = cumsum(Obs_inv,1);  