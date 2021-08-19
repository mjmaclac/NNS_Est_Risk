%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Import data figures and likelihood profiles for SI Appendix %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plotting data with missing info
temp1=imports_vol_id;
temp1(temp1==0)=NaN;
temp2=imports_vol_id;
temp2(temp1==1)=NaN;
temp2(isnan(temp1))=1;

yr_beg = 1940; 
yr_end = 1980;

figure(20); clf; 
figtitle = 'Fig_S1'; 
figh=2.5; figw=4.5; %fig width and height
set(gcf,'Units','inches','Position',[1,2,figw,figh]); %[left, bottom, width, height]
fs=11;  %fontsize

dataMsgInterpD = squeeze(temp2((yr_beg:yr_end)-1853,1,1)).*imports_vol_nan((yr_beg:yr_end)-1853,1,1);
dataReal       = squeeze(temp1((yr_beg:yr_end)-1853,1,1)).*imports_vol_nan((yr_beg:yr_end)-1853,1,1);
dataComb       = nansum([dataMsgInterpD dataReal],2);

plot(yr_beg:yr_end,dataComb,'b','LineWidth',1.5)
hold on
pl=plot(yr_beg:yr_end,dataMsgInterpD,'ro'); 
xlim([yr_beg yr_end])
set(gcf, 'Color', 'w')
xlabel('Year','FontSize',fs)
ylabel('Value of annual imports ($2015B)','FontSize',fs)
set(gca, 'XTick', [yr_beg:5:yr_end])
grid on
lg = legend(pl,'Missing/interpolated','Location','Northwest'); set(lg,'FontSize',fs)

print([export_dir figtitle],'-dpng','-r600')
     
%%    Value of imports -- annual and cumulative, aggregate (SI fig)
figure(30); clf; 
figtitle = 'Fig_S2'; 
figh=4.0; figw=6.5; %fig width and height
set(gcf,'Units','inches','Position',[1,5,figw,figh]); %[left, bottom, width, height]
fs=11;  %fontsize
xtickval = [1854 1900:50:2000];
    subplot(1,2,1)
        plot(1854:2012,squeeze(sum(Tot_val((1854:2012)-1853,use_regions,1),2))/1e9,'b','LineWidth',1.5)
        set(gca, 'XTick', xtickval )
        xlim([1854 2012])
        ylabel('Annual, \its_t', 'Fontsize',fs)
        grid on
        xlabel('Year, \itt', 'Fontsize',fs)
        text(1830,3.27,'A','Fontsize',12); 
    subplot(1,2,2)
        plot(1854:2012,cumsum(squeeze(sum(Tot_val((1854:2012)-1853,use_regions,1),2)))/1e9,'b','LineWidth',1.5)
        set(gcf, 'Color', 'w')
        xlabel('Year, \itt', 'Fontsize',fs)
        ylabel('Cumulative, \itS_t', 'Fontsize',fs)
        xlim([1854 2012])   
        grid on
        text(1825,75.5,'B','Fontsize',12); 
        set(gca, 'XTick', xtickval )
        h = axes('Position',[-.05 0 1 1],'Visible','off'); 
        set(gcf,'CurrentAxes',h)
        text(.08,.25,'Value of imports ($2015B)',...
            'VerticalAlignment','bottom',...
            'HorizontalAlignment','left', 'Rotation', 90, 'Fontsize',fs)        

print([export_dir figtitle],'-dpng','-r600')    

%% Negative log-likelihood profiles around the alpha parameters

figure(12)
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

print([export_dir 'Fig_S4'],'-dpng','-r600')