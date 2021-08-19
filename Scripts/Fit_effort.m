%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Native Hemiptera data and effort proxy %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

data=xlsread([spread_dir 'Data_out_basic']);                    % Two column count data
year=data(:,1)-min(data(:,1))+1;                                % Year is in the first column
count=log(cumsum(data(:,2)));                                   % Count is in the second column

N0=[5446.24946543026 0.0300533432587703 170.84690102395];       % Initial parameter guess
f2=@(N,year) log(N(1))-log((1+exp(-N(2).*(year-N(3)))));        % Three parameter logistic growth function
opts = statset('Display','iter','TolFun',1e-20,'TolX',1e-20);   % Options
z=fitnlm(year,count,f2,N0,'Options',opts);                      % Fitting program

t=1:(max(year));                                                % Time series for fitted logistic growth curve
count_est=f2(z.Coefficients.Estimate,t);                        % Fitted logistic growth curve

z.Coefficients

%% Make effort figure for SI (S.3)

figtitle = 'Fig_S3'; 

figure(1); clf
plot(t+1757,exp(count_est),'--r','LineWidth',1.5)               % Plot fitted curve
hold on 
plot(year+1757,exp(count),'b','LineWidth',1.5)                  % Plot data in red
xlim([1758 max(t+1757)]);grid on;                               % Reduce domain
set(gcf, 'Color', 'w')     
xlabel('Year')
ylabel({'Cumulative native','species discoveries'})
legend('Observed discoveries','Fitted discoveries', ...
    'Location','Northwest')

figh=6; figw=10;                                                % Fig width and height
set(gcf,'Units','inches','Position',[1,2,figw,figh]);       
set(findall(gcf,'-property','FontSize'),'FontSize',24)

print([export_dir figtitle],'-dpng','-r600')



%% Functional form for effort proxy

maxN = z.Coefficients.Estimate(1);                              % Max number of species
eofd1 = @(dt,Dt) -2500*log(1-dt./(maxN-Dt));                    % dt: discoveries; Dt: cumulative discoveries to t
eofd2 = @(dt,Dt) 2500*dt./(maxN-Dt); 
 
%% Smooth ITIS data: apportion discoveries in year t evenly over the specified window, "win"

win=10;
data_temp_ys =[1758-win+1:2012]';
data_temp =[data_temp_ys zeros(size(data_temp_ys))];
data_sm=data_temp;
for j=1:length(data(:,1))
   j_sm=find(data_temp(:,1)==data(j,1));
   data_temp(j_sm,2)=data(j,2);
end
for j=1:length(data_temp(:,1))
   if data_sm(j,1)>=2012-win+2
        win_temp=min(2012-data_sm(j,1)+1,win);
   else
        win_temp=win; 
   end
   data_sm(j,2) = sum((1/win_temp)*data_temp(j:(j+win_temp-1),2)); 
end

effort1=eofd1(data_sm(:,2),[0; cumsum(data_sm(1:end-1,2))]);
effort2=eofd2(data_sm(:,2),[0; cumsum(data_sm(1:end-1,2))]);

effort=effort2;
save([output_dir 'effort.mat'],'effort')