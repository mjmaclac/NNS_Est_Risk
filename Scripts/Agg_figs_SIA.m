%% %%%%%%%%%%
% Figure S5 %
%%%%%%%%%%%%%
clear text
temp_IMS   =cumsum(Invas_mean_est,1);
cum_IMS    = temp_IMS; % cumsum(Invas_mean_est,1);
cum_d      = cumsum(d_est,1);
Time_span2 = 1800:2012;

figure(13); clf; 
figtitle = 'Fig_S5'; 
figh=5; figw=6.5; %fig width and height
set(gcf,'Units','inches','Position',[1,2,figw,figh]); %[left, bottom, width, height]

fs       = 11;                                      % Fontsize
fsl      = 20;                                      % Large fontsize
kols     = 'krb';                                   % Three colors for three lines
symb     = ' . ';                                   % Symbols
lnsp     = {'-',':','--'};                          % Three line specs
ms       = 15;                                      % Marker size
mksp     = 20;                                      % Span between markers
xtickval = [1800 1854 1900:50:1950 2012];           % Displayed years on x-axis

% Plot estimated cumulative introductions versus discoveries
plot(Time_span2,sum(Cum_obs_inv(end,use_regions),2),[kols(1) lnsp{1} symb(1)],'MarkerSize',ms,'LineWidth',1.5,'MarkerIndices', 1:mksp:length(Time_span2)); hold on
plot(Time_span2,sum(cum_d,2)                      ,[kols(2) lnsp{2} symb(2)],'MarkerSize',ms,'LineWidth',1.5,'MarkerIndices', 1:mksp:length(Time_span2));
plot(Time_span2,sum(cum_IMS,2)                    ,[kols(3) lnsp{3} symb(3)],'MarkerSize',ms,'LineWidth',1.5,'MarkerIndices', 1:mksp:length(Time_span2));
xlim([1800 2016])
xlabel('Year')
ylabel('Number of species')
set(gca, 'XTick', []) 
set(gcf, 'Color', 'w')
set(findall(gcf,'-property','FontSize'),'FontSize',fs)
x0_0 = 2012;
x=2016;
y1 = sum(Cum_obs_inv(end,use_regions))+sum(Obs_inv_reg(end,use_regions));
y2 = sum(cum_IMS(end,:),2)+mean(Invas_mean_est(end,:));
width0=4;
width = 1;
line_x = x0_0  + [0, 0.5, 0.5, 1, 0.5, 0.5, 0]*width0;
line_y = y1 + [0, 0.02, 0.48, 0.5, 0.52, 0.98, 1]*(y2-y1);
line(line_x, line_y, 'Color', 'k')
text(x+width, y1 + 0.5*round(y2-y1), num2str(round(y2-y1)));
set(gca, 'XTick', xtickval)
box off
grid on
legend('Observed Discoveries','Fitted Discoveries','Estimated Establishments','location','northwest'); legend('boxoff')
set(findall(gcf,'-property','FontSize'),'FontSize',fs)

set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 0.4 0.7]);
         set(findall(gcf,'-property','FontSize'),'FontSize',fs)

figh=8; figw=10;                                        % Fig width and height
set(gcf,'Units','inches','Position',[1,2,figw,figh]);   % [left, bottom, width, height]
set(findall(gcf,'-property','FontSize'),'FontSize',18)
         
print([export_dir figtitle],'-dpng','-r600')

%% %%%%%%%%%%
% Figure S6 %
%%%%%%%%%%%%%

figtitle = 'Fig_S6'; 
reg_nums=1;
Invas_mean_fig= exp(parms(2)+repmat((parms(3)'),len,1).*log((exp(Tot_val_ln(:,1,3))+1e6)) ...
                            +repmat((parms(4)'),len,1).*Cum_Val_ln(:,1,3) ...
                            +parms(end)*repmat((-1:157)',1,length(reg_nums))); 
Time_span_temp=Time_span(2:end);
part_derivs_agg=Invas_mean_fig(2:end,:)-Invas_mean_est(Time_span_temp-1799,reg_nums);

% part_derivs_agg=repmat((parms(8+reg_nums)'),len,1).*Invas_mean_est(Time_span-1799,reg_nums)./squeeze(exp(Tot_val_ln(Time_span-1853,reg_nums,3)));
IME_hold=Invas_mean_est(end-(length(Time_span))+1:end,reg_nums);

sm_ln=10;
sm_cnt=floor(len/10)+1;

pd_yr=zeros(16,1);
IME_yr=zeros(16,1);
imp_yr=zeros(16,1);
for i=1:16
    j=10*i;
    if j==160
        j=158;
    end
    k=10*(i-1);
    hold1=cumsum(part_derivs_agg(1:j));
    hold3=cumsum(IME_hold(1:j));
    hold5=cumsum(Tot_val(1:j,1,3)*1e9);
    pd_yr(i,:)=hold1(end);
    IME_yr(i,:)=hold3(end);
    imp_yr(i,:)=hold5(end);
    if i==1
        pd_yr(i)=pd_yr(i)/sm_ln;
        IME_yr(i)=IME_yr(i)/sm_ln;
        imp_yr(i,:)=imp_yr(i)/sm_ln;
    elseif i<16
        hold2=cumsum(part_derivs_agg(1:k));
        pd_yr(i)=(pd_yr(i,:)-hold2(end))/sm_ln;
        hold4=cumsum(IME_hold(1:k));
        IME_yr(i)=(IME_yr(i)-hold4(end))/sm_ln;
        hold6=cumsum(Tot_val(1:k,1,3));
        imp_yr(i,:)=(imp_yr(i,:)-hold6(end,:))/sm_ln;
    else 
        hold2=cumsum(part_derivs_agg(1:k));
        pd_yr(i)=(pd_yr(i,:)-hold2(end))/(sm_ln-2);
        hold4=cumsum(IME_hold(1:k));
        IME_yr(i)=(IME_yr(i)-hold4(end))/(sm_ln-2);
        hold6=cumsum(Tot_val(1:k,1,3));
        imp_yr(i,:)=(imp_yr(i,:)-hold6(end,:))/(sm_ln-2);
    end
end
pd_yr=[pd_yr(1,:); pd_yr; pd_yr(end,:)];
IME_yr=[IME_yr(1,:); IME_yr; IME_yr(end,:)];
imp_yr=[imp_yr(1,:); imp_yr; imp_yr(end,:)];

ylabels={'Imports ($2015MM)','Marginal risk','Establishments'};
ts_yr=[1855 round(median(1855+sm_ln/2)) (1855+sm_ln):sm_ln:sm_ln*round((2012-sm_ln)/sm_ln)  round(median(2012-sm_ln/2)) 2012];

figure('Name','16');
        h=subplot(2,2,1);
        [AX,hlines]=plotyyy(h,ts_yr,imp_yr,ts_yr,pd_yr,ts_yr,IME_yr,ylabels); 
        xlim(AX(1), [min(ts_yr) max(ts_yr)]) 
        xlim(AX(2), [min(ts_yr) max(ts_yr)]) 
        tick_vec=[0 1 1e2 1e6 1e7 1e8 1e9];
        set (AX(1), 'Yscale', 'log','Ytick',tick_vec,'Yticklabels',string(tick_vec*1e-6));
        set (AX(1), 'Xtick',[1855 1900 1950 2012],'Xticklabels',string([1855 1900 1950 2012]));
        xlabel(AX(1),'Year')
        set(gcf, 'Color', 'w')   
        set(AX(3),'xcolor','none')
        hL = legend(hlines,'Orientation','horizontal','Imports (left axis)','Marginal risk (right axis)','Establishments (right axis)'); %,'Location', 'southoutside'
        legend boxoff
        newPosition = [0.25 0.4 0.1 0.1];
        newUnits = 'normalized';
        set(hL,'Position', newPosition,'Units', newUnits);        
        
figh=12; figw=24;                                      % Fig width and height
set(gcf,'Units','inches','Position',[1,2,figw,figh]);  % [left, bottom, width, height]
set(findall(gcf,'-property','FontSize'),'FontSize',24)

print([export_dir figtitle],'-dpng','-r600')

%% %%%%%%%%%%
% Figure S7 %
%%%%%%%%%%%%%

figtitle = 'Fig_S7'; 

figure(17)
    subplot(3,1,1)
        plot(a_0,-squeeze(min(min(fval_array_temp,[],3),[],2)))
        title('\alpha_0', 'Interpreter', 'tex')
        box off
        grid on
    subplot(3,1,2)
        plot(a_1,-squeeze(min(min(fval_array_temp,[],3),[],1)))
        title('\alpha_1', 'Interpreter', 'tex')
        box off
        grid on
    subplot(3,1,3)
        plot(a_2,-squeeze(min(min(fval_array_temp,[],2),[],1)))
        title('\alpha_2', 'Interpreter', 'tex')
        box off
        grid on
        set(gcf, 'Color', 'w')
 
print([export_dir figtitle],'-dpng','-r600')