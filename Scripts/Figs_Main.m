%% %%%%%%%%%
% Figure 3 %
%%%%%%%%%%%%
clear text
temp_IMS   = cumsum(Invas_mean_est,1);
cum_IMS    = temp_IMS; % cumsum(Invas_mean_est,1);
cum_d      = cumsum(d_est,1);
Time_span2 =1800:2012;

fig        = figure(1); clf; 
figtitle = 'Fig_3'; 
figh=10; figw=20;       % Fig width and height

fs   = 30;              % Fontsize basic
fsl  = 30;              % Fontsize large (same as above)
kols = 'krb';           % Three colors for three lines
symb = ' . ';           % Three symbols
lnsp = {'-',':','--'};  % Three line specs
ms   = 15;              % Marker size
mksp = 20;              % Span between markers
xtickval = [1800 1855 1900:50:1950 2012];

make_it_tight=1;
subplot = @(m,n,p) subtightplot (m, n, p, [0.12 0.175], [0.125 0.1], [0.1 0.025]);
if ~make_it_tight,  clear subplot;  end

% Plot estimated cumulative introductions versus discoveries
subplot (1,2,1)
    box off
    p1 = plot(Time_span2,sum(Cum_obs_inv(:,use_regions),2)',[kols(1) lnsp{1} symb(1)],'MarkerSize',ms,'LineWidth',1.5,'MarkerIndices', 1:mksp:length(Time_span2)); hold on
    p2 = plot(Time_span2,sum(cum_d,2)                      ,[kols(2) lnsp{2} symb(2)],'MarkerSize',ms,'LineWidth',1.5,'MarkerIndices', 1:mksp:length(Time_span2));
    p3 = plot(Time_span2,sum(cum_IMS,2)                    ,[kols(3) lnsp{3} symb(3)],'MarkerSize',ms,'LineWidth',1.5,'MarkerIndices', 1:mksp:length(Time_span2));
    xlim([1800 2018])
    xlabel('Year')
    ylabel('Number of species')
    set(gcf, 'Color', 'w')
    set(gca, 'XTick', xtickval)
    x0_0 = 2012;
    x=2016;
    y1 = sum(Cum_obs_inv(end,use_regions));
    y2 = sum(cum_IMS(end,:),2)+sum(Invas_mean_est(end,:));
    width0= 6;
    width = 4;
    line_x = x0_0  + [0, 0.5, 0.5, 1, 0.5, 0.5, 0]*width0;
    line_y = y1 + [0, 0.02, 0.48, 0.5, 0.52, 0.98, 1]*(y2-y1);
    line(line_x, line_y, 'Color', [0.4940 0.1840 0.5560],'LineWidth',6)
    text(x+width, y1 + 0.5*round(y2-y1), {num2str(round(y2-y1));join(['(' string(round(100*(y2-y1)/y2)) '%)'],"")},'Color',[0.4940 0.1840 0.5560]);
    line_y = y1 + [0, 0.02, 0.48, 0.5, 0.52, 0.98, 1]*(-y1);
    text(x+width, y1 + 0.5*round(0-y1), num2str(round(y1)),'Color',[0 .5 .5]);
    line(line_x, line_y, 'Color', [0 .5 .5],'LineWidth',6)
    box off
    set(gca,'TickLength',[0 0])
    grid on
    legend([p1 p2 p3],{'Observed Discoveries','Fitted Discoveries','Estimated Establishments'},'location','northwest'); legend('boxoff')
    ylimtemp = get(gca,'YLim');
    text(1800,1.03*ylimtemp(2),'(A)','HorizontalAlignment','right','VerticalAlignment','bottom','FontSize',fsl,'fontweight','bold');
    text(2125,1.03*ylimtemp(2),'(B)','HorizontalAlignment','right','VerticalAlignment','bottom','FontSize',fsl,'fontweight','bold');
 subplot (1,2,2)
    X=categorical({'European';'Neotropic';'Asian';'Australasia';'Indomalaya';'Afrotropic';'Oceania'});
    X=reordercats(X,{'European';'Neotropic';'Asian';'Australasia';'Indomalaya';'Afrotropic';'Oceania'});
    Y=[Cum_obs_inv(end,use_regions([4 6 2 3 5 1 7]));cum_IMS(end,[4 6 2 3 5 1 7])-Cum_obs_inv(end,use_regions([4 6 2 3 5 1 7]))]';
        bb=barh(X,Y,'stacked','FaceColor',[0 .5 .5],'EdgeColor',[0 .9 .9],'LineWidth',0.5);
        bb(2).FaceColor=[0.4940 0.1840 0.5560];
        bb(2).FaceAlpha=0.4;
        xlabel('Number of species')
        legend(bb,'Discovered species','Establishment debt')  
        legend('boxoff')
        labelArray = {'European Palearctic';'Neotropic';'Asian Palearctic';'Australasia';'Indomalaya';'Afrotropic';'Oceania'};
        labels = cellfun(@(x) strrep(x,' ','\newline'), labelArray ,'UniformOutput',false);
        set(gca,'TickLength',[0 0])
        yticklabels(labels)
        ax=gca;
        ax.XGrid = 'on';
        ax.YGrid = 'off';
        yt = get(gca, 'YTick');
        labels=join([string(round(100*(cum_IMS(end,[4 6 2 3 5 1 7])-Cum_obs_inv(end,use_regions([4 6 2 3 5 1 7])))./ ...
        cum_IMS(end,[4 6 2 3 5 1 7])))',repmat('%',7,1)],"");
        x_dist=cum_IMS(end,[4 6 2 3 5 1 7]);
        for i=1:7
            if i>2
               text(x_dist(i)+10, yt(i),labels(i), 'HorizontalAlignment','left','Color',[0.4940 0.1840 0.5560])    
            else
               text(x_dist(i)-10, yt(i),labels(i), 'HorizontalAlignment','right','Color',[0.4940 0.1840 0.5560]) 
            end
        end
   
    box off
    set(findall(gcf,'-property','FontSize'),'FontSize',fs)
    set(gcf,'Units','inches','Position',[1,2,figw,figh]); %[left, bottom, width, height]
             
print([export_dir figtitle],'-dpng','-r600')

%% %%%%%%%%%
% Figure 4 %
%%%%%%%%%%%%

figure(4)
reg_nums=1:7;
Invas_mean_fig= exp(parms(8)+repmat((parms(8+reg_nums)'),len,1).*log((exp(Tot_val_ln(:,use_regions(reg_nums),3))+1e6)) ...
                            +repmat((parms(15+reg_nums)'),len,1).*Cum_Val_ln(:,use_regions(reg_nums),3) ...
                            +parms(end)*repmat((-1:157)',1,length(reg_nums))); 
Time_span_temp=Time_span(2:end);
part_derivs_agg=Invas_mean_fig(2:end,:)-Invas_mean_est(Time_span_temp-1799,reg_nums);
  
IME_hold=Invas_mean_est(end-(length(Time_span))+1:end,reg_nums);

sm_ln=10;
sm_cnt=floor(len/10)+1;

pd_yr=zeros(sm_cnt,7);
IME_yr=zeros(sm_cnt,7);
imp_yr=zeros(sm_cnt,7);
for i=1:sm_cnt
    j=sm_ln*i;
    if j==160
        j=158;
    end
    k=sm_ln*(i-1);
    hold1=cumsum(part_derivs_agg(1:j,:));
    hold3=cumsum(IME_hold(1:j,:));
    hold5=cumsum(Tot_val(1:j,use_regions(reg_nums),3));
    pd_yr(i,:)=hold1(end,:);
    IME_yr(i,:)=hold3(end,:);
    imp_yr(i,:)=hold5(end,:);
    if i==1
        pd_yr(i,:)=pd_yr(i,:)/sm_ln;
        IME_yr(i,:)=IME_yr(i,:)/sm_ln;
        imp_yr(i,:)=imp_yr(i,:)/sm_ln;
    elseif i<sm_cnt 
        hold2=cumsum(part_derivs_agg(1:k,:));
        pd_yr(i,:)=(pd_yr(i,:)-hold2(end,:))/sm_ln;
        hold4=cumsum(IME_hold(1:k,:));
        IME_yr(i,:)=(IME_yr(i,:)-hold4(end,:))/sm_ln;
        hold6=cumsum(Tot_val(1:k,use_regions(reg_nums),3));
        imp_yr(i,:)=(imp_yr(i,:)-hold6(end,:))/sm_ln;
    else
        hold2=cumsum(part_derivs_agg(1:k,:));
        pd_yr(i,:)=(pd_yr(i,:)-hold2(end,:))/(sm_ln-2);
        hold4=cumsum(IME_hold(1:k,:));
        IME_yr(i,:)=(IME_yr(i,:)-hold4(end,:))/(sm_ln-2);
        hold6=cumsum(Tot_val(1:k,use_regions(reg_nums),3));
        imp_yr(i,:)=(imp_yr(i,:)-hold6(end,:))/(sm_ln-2);
    end
end
pd_yr=[pd_yr(1,:); pd_yr; pd_yr(end,:)];
IME_yr=[IME_yr(1,:); IME_yr; IME_yr(end,:)];
imp_yr=[imp_yr(1,:); imp_yr; imp_yr(end,:)];

make_it_tight = true;
subplot = @(m,n,p) subtightplot (m, n, p, [0.075 0.13], [0.1 0.05], [0.1 0.1]);
if ~make_it_tight,  clear subplot;  end
figtitle='Fig_4';

graph_ord=[6 4 2 5 3 1];
ts_yr=[1855 round(median(1855+sm_ln/2)) (1855+sm_ln):sm_ln:sm_ln*round((2012-sm_ln)/sm_ln)  round(median(2012-sm_ln/2)) 2012];
letr = 'ABCDEF';

fig=figure('Name','5','units','normalized','outerposition',[0 0 0.4 0.6]);
    for i=1:6
            h=subplot(4,2,i);
            j=graph_ord(i);
            if j>5
                j=j+1;
            end
            ylabels={'Imports ($2015MM)','Marginal risk','Establishments'};
            [AX,hlines]=plotyyy(h,ts_yr,imp_yr(:,graph_ord(i)),ts_yr,pd_yr(:,graph_ord(i)),ts_yr,IME_yr(:,graph_ord(i)),ylabels); 
            xlim(AX(1), [min(ts_yr) max(ts_yr)]) 
            xlim(AX(2), [min(ts_yr) max(ts_yr)]) 
            set(AX(3),'xcolor','none')
            title(reg_list(j+1))
            if i<5
                set(AX, 'XTick', []) 
            else
                set(AX, 'XTick', [min(ts_yr) 1900 1950 max(ts_yr)])  
                xlabel(AX(1),'Year')
            end
            if i==1 || i==3 || i==5
               ylabel(AX(2),'')
               ylabel(AX(3),'')
            else
               ylabel(AX(1),'')
            end
            set(gcf, 'Color', 'w') 
            ylimtemp = get(gca,'YLim');
            text(1855,1.05*ylimtemp(2),['(' [letr(i)] ')'],'HorizontalAlignment','right','VerticalAlignment','bottom','FontSize',fsl,'fontweight','bold');
            if i==5 || i==6; xlabel('Year','FontSize',fsl); end
            if i==6
                hL = legend(hlines,'Orientation','horizontal',{'Imports (left axis)','Marginal risk (right axis)','Establishments (right axis)'});
                legend boxoff
                newPosition = [0.5 0.16 0.1 0.1];
                newUnits = 'normalized';
                set(hL,'Position', newPosition,'Units', newUnits);
            end 
    end       
set(gcf,'Units','inches','Position',[1,2,figw,1.33*figh]); %[left, bottom, width, height]
set(findall(gcf,'-property','FontSize'),'FontSize',20)

print([export_dir figtitle],'-dpng','-r600')

%% %%%%%%%%%
% Figure 5 %
%%%%%%%%%%%%

close all
subplot = @(m,n,p) subtightplot (m, n, p, [0.15 0.13], [0.1 0.1], [0.1 0.1]);
if ~make_it_tight,  clear subplot;  end
figtitle='Fig_5';
letr = 'AB';

figure(5)
    subplot(2,1,1)
        plot(Time_span2,diag(put),'LineWidth',2)
        xlim([1800 2012])
        ylabel('Initial probability of discovery')
        xlabel('Year')
        ylimtemp = get(gca,'YLim');
        text(1800,1.05*ylimtemp(2),['(' [letr(1)] ')'],'HorizontalAlignment','right','VerticalAlignment','bottom','FontSize',fsl,'fontweight','bold');
    subplot(2,1,2)
        plot(0:100,put(1:101,1),'-','LineWidth',2); hold on
        plot(0:100,put(51:151,51),'--','LineWidth',2); hold on
        plot(0:100,put(101:201,101),':','LineWidth',2); hold on
        plot(0:62,put(151:213,151),'-*','LineWidth',1); 
        ylabel('Probability of discovery')
        xlabel('Years after establishment')
        legend('1800','1850','1900','1950','location','Northwest')
        set(gcf, 'Color', 'w')
        ylimtemp = get(gca,'YLim');
        text(0,1.2*ylimtemp(2),['(' [letr(2)] ')'],'HorizontalAlignment','right','VerticalAlignment','bottom','FontSize',fsl,'fontweight','bold');

figh=12; figw=20; %fig width and height
set(gcf,'Units','inches','Position',[1,2,figw,figh]); %[left, bottom, width, height]
set(findall(gcf,'-property','FontSize'),'FontSize',24)

print([export_dir figtitle],'-dpng','-r600')