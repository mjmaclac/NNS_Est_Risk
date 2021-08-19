%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Model of trade, introduction and discovery %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Early import data (1854-1966) (do not run if import_data is loaded)
[trash,conve]=xlsread([spread_dir 'GATS_Basic'],'Partner_Region'); clear trash;
reg_list={'Aggregate','Afrotropic','Asian Palearctic','Australasia','European Palearctic','Indomalaya','Nearctic','Neotropic','Oceania','Misc'};
reg_list_b={'Aggregate','Afrotropic','Asian P.',         'Austral.','Euro. P.',         'Indomal.',      'Nearctic','Neotropic','Oceania','Misc'};
Prod_list=['A','S','N','C'];
% A="Aggregate"; S="Seeds"; N="Nursery products"; C="Cut flowers"
imports_vol=zeros(2014-1854+1,9,4);  % Initialize array of later import values 

for t=[1854:1946 1949:1965]   % Data are available for 1854 - 1965 (with gaps); only seeds are available for 1966. 
    clear eco_reg;
    if (t>1949 && t<1956)  || t==1957
       M=1:12;
    else
       M=1; 
    end
    for mm=M
        if t<1931
            [num, txt]=xlsread([spread_dir 'Trade_data_early.xlsx'], num2str(t));
        elseif (t>1949 && t<1956)  || t==1957
            [num, txt]=xlsread([spread_dir 'Trade_data_late.xlsx'], [num2str(t) '_' num2str(mm)]);
        else
            [num, txt]=xlsread([spread_dir 'Trade_data_late.xlsx'], num2str(t));
        end
        num(isnan(num))=0;
        eco_reg(1,1)={'Aggregate'};
        for i=3:size(txt,1)
           for j=1:size(conve,1) 
              if strcmp(txt(i,1),conve(j,1))
                  eco_reg(i-1,1)=conve(j,2);
              end
           end
        end
        for i=1:size(num,1)
            for r=1:9
                for p=1:4
                    for m=2:size(txt,2)
                        if strcmp(eco_reg(i,1),reg_list(r)) && strncmpi(Prod_list(p),txt(1,m),1) 
                            imports_vol(t-1853,r,p)=max(imports_vol(t-1853,r,p),0)+max(num(i,m-1),0);
                        end
                    end
                end
            end
        end
    end
end

% Late import data (1967-2012)
[num, txt]=xlsread([spread_dir 'Data_by_prod.xls'], 'Sheet1');
num(:,2)=num(:,2)*1e6;
for t=1967:2012
    for i=1:size(num,1)
        for r=1:9
            for p=1:4
                if t==num(i,1) && strncmpi(Prod_list(p),txt(i,2),1) && strcmp(txt(i,1),reg_list(r))
                    imports_vol(t-1853,r,p)=max(imports_vol(t-1853,r,p),0)+num(i,2);
                end
            end
        end
    end
end

%% Summing across product groups 
imports_vol((1967:2012)-1853,1,:)=sum(imports_vol((1967:2012)-1853,:,:),2);
imports_vol(:,:,1)=sum(imports_vol,3);
imports_vol=imports_vol/1e9;

%% Adjust import value data for inflation
[num, txt]=xlsread([spread_dir 'inflation.xlsx'], 'Sheet4');
for t=(1854:2012)-1853
    for r=1:9
        for p=1:4
            imports_vol(t,r,p)=imports_vol(t,r,p)*num(t,2);
        end
    end
end
imports_vol_id=ones(size(imports_vol));
imports_vol_uninterp=imports_vol;

% Aggregate interpolation 
for t=1854:1966
    skip=0;
    if imports_vol(t-1853,1,1)==0
        while imports_vol(t-1853+skip,1,1)==0 && skip+t<2013
            skip=skip+1;
        end
        imports_vol(t-1853,1,1)=((skip)/(1+skip))*imports_vol(t-1853-1,1,1)+(1/(1+skip))*imports_vol(t-1853+skip,1,1);
        imports_vol_id(t-1853,1,1)=NaN;
    end
end

%% Regional interpolation 
reg_prop=imports_vol(:,:,1)./repmat(imports_vol(:,1,1),[1,9,1]);
reg_prop(isnan(reg_prop))=0;

for t=1854:1966
    for r=2:9    
        skip=1;
        if imports_vol(t-1853,1,1)*0.9<=sum(imports_vol(t-1853,2:9,1),2) && round(1e9*imports_vol(t-1853,1,1))>round(1e9*sum(imports_vol(t-1853,2:9,1),2)) 
            imports_vol(t-1853,r,1)=reg_prop(t-1853,r).*imports_vol(t-1853,1,1)./sum(reg_prop(t-1853,:));
        elseif imports_vol(t-1853,1,1)*0.9>sum(imports_vol(t-1853,2:9,1),2) 
            while sum(reg_prop(t-1853+skip,2:9),2)<1 && skip+t<2013
                skip=skip+1;
                reg_prop(t-1853,r)=(((skip)/(1+skip)).*reg_prop(t-1853-1,r)+(1/(1+skip)).*reg_prop(t-1853+skip,r));
            end 
        end
    end
    reg_prop(t-1853,2:9)=reg_prop(t-1853,2:9)./sum(reg_prop(t-1853,2:9));
    imports_vol(t-1853,2:9,1)=reg_prop(t-1853,2:9)*imports_vol(t-1853,1,1);
end
    imports_vol(isnan(imports_vol))=0;

%% Interpolating missing product data
prod_prop=imports_vol(:,:,:)./repmat(imports_vol(:,:,1),[1,1,4]);
prod_prop(isnan(prod_prop))=0;

for t=1854:1966
    for r=1:9   
        for p=2:4
            skip=1;
            if t>1854 && imports_vol(t-1853,r,1)>nansum(imports_vol(t-1853,r,2:4),3) 
                while prod_prop(t-1853+skip,r,p)==0 && skip+t<2013
                   skip=skip+1;
                end
                prod_prop(t-1853,r,p)=(((skip)/(1+skip)).*prod_prop(t-1853-1,r,p)+(1/(1+skip)).*prod_prop(t-1853+skip,r,p));
            elseif imports_vol(t-1853,r,1)>nansum(imports_vol(t-1853,r,2:4),3) 
                while prod_prop(t-1853+skip,r,p)==0 && skip+t<2013
                   skip=skip+1;
                end
                prod_prop(t-1853,r,p)=prod_prop(t-1853+skip,r,p);
             end
        end
        if imports_vol(t-1853,r,1)>0
            prod_prop(t-1853,r,2:4)=prod_prop(t-1853,r,2:4)./sum(prod_prop(t-1853,r,2:4));
            imports_vol(t-1853,r,2:4)=prod_prop(t-1853,r,2:4)*imports_vol(t-1853,r,1);
        else 
            imports_vol(t-1853,r,2:4)=zeros(3,1);
        end
    end
end

%% Modifications
imports_vol_nan=imports_vol(1:end-2,:,:);       % Data end in year 2012
imports_vol_nan(isnan(imports_vol_nan))=0;      % Converts any NaN entries to 0 

% Defining total value (Tot_val) variable
Tot_val=squeeze(imports_vol_nan(:,:,1));
Tot_val(isnan(Tot_val))=0;

%%
Cum_Val_hold=cumsum(Tot_val);
Cum_Val=[zeros(1,9); Cum_Val_hold(1:end-1,:)];    % Calculate cumulative value  