%Calculate the WPE of the logistic map across various levels of Tau and Permutation length

clearvars
close all
cd('/Users/alisoniles/Documents/GitHub/Data/Logistic-Map/Logistic_map_repo');

M = 10000; % number of iterations of logistic equation
N = 3000; % length of end of time series to save and analyzes

R = 500; % number of "r" values to simulate
a = 3.4; % starting value of "r"
b = 3.9; % final value of "r"... anything higher diverges.
    rs = linspace(a,b,R); % vector of "r" values

mwl=5; % Maximum word length of permutations to consider
mtau=4; % Maximum time delay between the 'letters' of each 'word' to consider

scale = 10000; % determines the level of rounding for plotting the bifurcation diagram (not for analysis)
maxpoints = 200; % determines maximum values to plot

pn = [0;0.0001;0.001;0.01]; % standard deviation of the process noise (added Gaussian dynamical noise to the state variable, x, at each time step) levels to loop through
WPE_Zpn=[]; WPE_Lpn=[]; WPE_Mpn=[]; WPE_Hpn=[];% allocate matrices to store WPE output for zero, low, medium and high process noise
LM_Zpn=zeros(N,R); LM_Lpn=zeros(N,R); LM_Mpn=zeros(N,R); LM_Hpn=zeros(N,R); % allocate matrices to store LM time series output for zero, low, medium and high process noise
xdiff_Zpn=zeros(M,R); xdiff_Lpn=zeros(M,R); xdiff_Mpn=zeros(M,R); xdiff_Hpn=zeros(M,R); % allocate matrices to store LM time series output for zero, low, medium and high process noise

%Loop throught the levels of process noise 'pn'
for k=1:length(pn)
    
% Loop through the "r" values
for j = 1:length(rs)
    
    r=rs(j); % get current "r"
    x=zeros(M,1); % allocate memory
    xdiff=zeros(M,1);
    x(1) = 0.5; % initial condition (can be anything from 0 to 1)
    xdiff(1) = 0; 
    
    for i = 2:M % iterate the logistic map with process noise
        
        %Add process noise to x at every time step, analogous to migration
        %becasuse the growth rate is not changing, just the population
        %size.   
        
        x(i) = r*x(i-1)*(1-x(i-1))+(pn(k).*randn(1,1)); %Logistic map with error. if pn>0, adds process noise to the x at every time step
          if x(i)<=0 || x(i)>=1 %If the new x value is not between 0 and 1, then
              NEWxList=r*x(i-1)*(1-x(i-1))+pn(k).*randn(100,1); %Generate a list of 100 NEWx values to chose the first one between 0 and 1
            NEWx=NEWxList(NEWxList>0 & NEWxList<1);
            x(i)=NEWx(1,1);
          end
         
         %summarize the resulting noise distribution
         x_expct=r*x(i-1)*(1-x(i-1)); %the expected x with no noise
         xdiff(i)=x(i)-x_expct; %the difference between x and the expected x
        
    end
    
    % Calculate the Weighted Permutation Entropy for process noise
        % Loop through the values of tau
        x(1:M-N,:)=[]; %Remove initial transient dynamics
        
        xr=round(x, 10); %Round to nearest 10^-10 to remove minor rounding fluctuations
        for t=1:mtau
            tau=t;
        [ H ] = WPE( xr, mwl, tau );
        
            if k==1; WPE_Zpn=[WPE_Zpn; r, tau, H, pn(k)];
            elseif k==2; WPE_Lpn=[WPE_Lpn; r, tau, H, pn(k)];
            elseif k==3; WPE_Mpn=[WPE_Mpn; r, tau, H, pn(k)];
            elseif k==4; WPE_Hpn=[WPE_Hpn; r, tau, H, pn(k)];
            end
        end
    
    % only save those unique, semi-stable values from zero process
    % noise to plot in bifurcation plot.
        if k==1; out{j} = unique(round(scale*x(end-maxpoints:end))); end % only save bifurcation data for zero added noise otherwise the plot is messy and unrecognizable as the logistic map
    
    % save LM time series data for export to R
        if k==1; LM_Zpn(:,j)=xr; xdiff_Zpn(:,j)=xdiff;
        elseif k==2; LM_Lpn(:,j)=xr; xdiff_Lpn(:,j)=xdiff;
        elseif k==3; LM_Mpn(:,j)=xr; xdiff_Mpn(:,j)=xdiff;
        elseif k==4; LM_Hpn(:,j)=xr; xdiff_Hpn(:,j)=xdiff;
        end
        
        
end %r values loop

end %process noise loop

% Adding observational noise
    on = [0; 0.0001; 0.001; 0.01]; % standard deviation of the observational noise (added Gaussian dynamical noise to the time series of x at the end of the simulation) levels to loop through

    LM_Zon=(on(1).*LM_Zpn).*randn(N,R)+LM_Zpn;
    LM_Lon=(on(2).*LM_Zpn).*randn(N,R)+LM_Zpn;
    LM_Mon=(on(3).*LM_Zpn).*randn(N,R)+LM_Zpn;
    LM_Hon=(on(4).*LM_Zpn).*randn(N,R)+LM_Zpn;
    
    % Calculate the Weighted Permutation Entropy for observational noise
      WPE_Zon=[]; WPE_Lon=[]; WPE_Mon=[]; WPE_Hon=[];  % allocate matrices to store WPE output for zero, low, medium, high and superhigh observational noise

            for j = 1:length(rs) % Loop through the different r values
                r=rs(j); 
                
            for t=1:mtau % Loop through the values of tau for WPE calculation
                tau=t;    
                   [ H ] = WPE(LM_Zon(:,j), mwl, tau);
                   WPE_Zon=[WPE_Zon; r, tau, H, on(1)];

                   [ H ] = WPE(LM_Lon(:,j), mwl, tau);
                   WPE_Lon=[WPE_Lon; r, tau, H, on(2)];

                   [ H ] = WPE(LM_Mon(:,j), mwl, tau);
                   WPE_Mon=[WPE_Mon; r, tau, H, on(3)];

                   [ H ] = WPE(LM_Hon(:,j), mwl, tau);
                   WPE_Hon=[WPE_Hon; r, tau, H, on(4)];

            end %tau loop
            end %r loop    
        
save 'LogisticMapTimeSeriesData.mat' LM_Zpn LM_Lpn LM_Mpn LM_Hpn LM_Zon LM_Lon LM_Mon LM_Hon rs; %Save time series so to access from R and calcuate the S-map projection and the nRMSE
save 'LogisticMapWPEData.mat' WPE_Zpn WPE_Lpn WPE_Mpn WPE_Hpn WPE_Zon WPE_Lon WPE_Mon WPE_Hon; %Save WPE data


%% Figure WPE vs r

% Rearrange cell array into a large n-by-2 vector for plotting
bfrctn = [];
for k = 1:length(rs)
    n = length(out{k});
    bfrctn = [bfrctn;  rs(k)*ones(n,1),out{k}];
end

% Plot the data
WPE_data=WPE_Zpn;% plot the zero process noise results
close all
f=figure(1);clf
hold all
set(f, 'Position', [100, 100, 1000, 1000]);
set(f, 'DefaultAxesFontSize', 12);
set(f, 'DefaultLineLineWidth',1.2);
set(f, 'color','white')
Alphabet = 'abcdefghijkl';

for t=1:mtau
    for w=3:mwl
        p=((t-1)*3)+w-2;
        sp=subplot(4,3,p);
        
        h1=plot(bfrctn(:,1),bfrctn(:,2)/scale,'k.');
        set(h1,'markersize',1)
        hold on

        h2=plot(WPE_data(WPE_data(:,2)==t,1),WPE_data(WPE_data(:,2)==t,w),'ro');
        set(h2,'markersize',5)
        hold on
        
        text(2.1,0.85,strcat(Alphabet(p),') \tau=',num2str(t),', wl=', num2str(w)));
        
        set(gca, ...           
            'XTick', 2:0.5:4           , ...
            'XLim', [a-0.1, b+0.1]         , ...
            'YTick', 0:0.2:1      , ...
            'xticklabel',[] , ...
            'yticklabel',[] , ...
            'YLim', [0, 1]    , ...
            'Box', 'on'             , ...
            'Linewidth', 1          , ...
            'TickDir','out')

        if p==1 || p==4 || p==7 || p==10, set(gca, 'yticklabel',{'0.0','0.2','0.4','0.6','0.8','1.0'}),end
        if p==10 || p==11 || p==12 , set(gca, 'xticklabel',{'2.0','2.5','3.0','3.5','4.0'}),end

        %Modify the positions/size of each subplot so that there is less
        %white space
        pos = get(sp, 'position'); %Get positions of each subplot
        if p==2 || p==5 || p==8 || p==11,  pos(1)=pos(1)+0.03; end
        if p==1 || p==4 || p==7 || p==10,   pos(1)=pos(1)+0.06; end
        if p==4 || p==5 || p==6,   pos(2)=pos(2)+0.03;  end
        if p==7 || p==8 || p==9,    pos(2)=pos(2)+0.06; end
        if p==10 || p==11 || p==12,   pos(2)=pos(2)+0.09;  end
        set(sp, 'position', pos);
        
        if p==11, xlabel('Growth rate, r'); end
        if p==7 
            yL=ylabel('WPE'); 
            pyL = get(yL, 'position'); %Modify position of y-axis label
            set(yL, 'position', [pyL(1),(pyL(2)+0.5),pyL(3)])
        end
    end
end

%% Get forecast error

% This part of the analysis is done in R. Load into R the Matlab data file 
% of the logistic map time series ('LogisticMapTimeSeriesData.mat') using the R script
% 'EDM_LogisticMap_noise_for_Matlab_v2.R, which uses the function 'nRMSE.R' to calculate the
% nRMSE of the logistic map timeseries at each of 100 
% projected time steps using the previous 1900 timestep to formulate the 
% attractor of the S-Map. 
 

%% Main Logistic Map Summary figure (figure 4 in MS)

load LogisticMapForecastError.mat; %load forecast error results from R back into MatLab
RMSE=output.rmse(strcmp(output.noise_level,'0_zero') & strcmp(output.noise_type,'process_noise'));

f4=figure(4);clf
hold all
set(f4, 'Position', [100, 100, 400, 800]);
set(f4, 'DefaultAxesFontSize', 14);
set(f4, 'DefaultLineLineWidth',1.5);
set(f4, 'color','white')

t=1; %The chosen tau to plot

subplot(4,1,1) %Demonstrate the conceptual leap from time series to bifurcation diagram. Plot 5 different time series of 3 different growth rates in the corresponding colors.

% Choose growth rate values to highlight in figures
r_plot=[3.4140280561 3.7266533066 3.8398797595];

col=colormap;

xvals=linspace(2970,3000,30)';

arcol=zeros(length(r_plot),3);

for l=1:length(r_plot)
    ll=find(round(rs,10)==r_plot(l),1); 
    yvals=LM_Zpn(2971:3000,ll);
    curcol=round(((length(col)-1)/length(rs))*ll)+1; %point color on plot
    arcol(l,:)=col(curcol,:);
    h1=plot(xvals,yvals, 'o-', 'Color',col(curcol,:), 'Linewidth', 1.5, 'MarkerFaceColor',col(curcol,:));
        hold on   
end

 set(gca, ... 
            'XLim', [N-30, N]         , ...
            'YLim', [-0.05, 1.05]    , ...
            'YTick', 0:0.2:1      , ...
            'yticklabel',{'0','0.2','0.4','0.6','0.8','1'} , ...
            'Linewidth', 1, ...
            'Box','on')
        
        xlabel('time step')
        ylabel('Abundance'); 
        text(2965.5,1.2,'A)','FontSize', 12, 'FontWeight', 'Bold');


subplot(4,1,2)
        h3=scatter(bfrctn(1:end-1,1),bfrctn(1:end-1,2)/scale,2,bfrctn(1:end-1,1),'o','filled');
        hold on
       
        set(gca, ... 
            'XLim', [a, b]         , ...
            'XTick', 3.4:0.1:3.9      , ...
            'xticklabel',{'3.4','3.5','3.6','3.7','3.8','3.9'} , ...
            'YLim', [-0.05, 1.05]    , ...
            'YTick', 0:0.2:1      , ...
            'yticklabel',{'0','0.2','0.4','0.6','0.8','1'} , ...
            'Linewidth', 1, ...
            'Box','on')
        
         for ii=1:length(r_plot)
            plot([r_plot(ii),r_plot(ii)],[1.05,-0.05], '.-', 'Color',arcol(ii,:), 'Linewidth', 2);
        hold on   
         end 
         
        xlabel('Growth rate, r')
        ylabel('Abundance'); 
        text(3.325,1.2,'B)','FontSize', 12, 'FontWeight', 'Bold');

subplot(4,1,3) 
WPE_data=WPE_Zpn;% plot the zero process noise results

plot(WPE_data(WPE_data(1:end-4,2)==t,1),WPE_data(WPE_data(1:end-4,2)==t,3),'-', 'Color',[0.7,0.7,0.7],'Linewidth', 1.5); hold on
plot(WPE_data(WPE_data(1:end-4,2)==t,1),WPE_data(WPE_data(1:end-4,2)==t,4),'-', 'Color',[0.4,0.4,0.4],'Linewidth', 1.5); hold on
plot(WPE_data(WPE_data(1:end-4,2)==t,1),WPE_data(WPE_data(1:end-4,2)==t,5),'k-', 'Linewidth', 1.5); hold on

      set(gca, ... 
            'XLim', [a, b]         , ...
            'XTick', 3.4:0.1:3.9      , ...
            'xticklabel',{'3.4','3.5','3.6','3.7','3.8','3.9'} , ...
            'YLim', [-0.05, 1.05]    , ...
            'YTick', 0:0.2:1      , ...
            'yticklabel',{'0','0.2','0.4','0.6','0.8','1'} , ...
            'Linewidth', 1, ...
            'Box','on')
        
            
        for ii=1:length(r_plot)
            iiy=find(round(rs,10)==r_plot(ii),1); 
            iix=find(round(WPE_data(:,1),10)==r_plot(ii) & WPE_data(:,2)==t);
            scatter(r_plot(ii),WPE_data(iix,5), 60, arcol(ii,:),'o','filled')
        end
        
        xlabel('Growth rate, r')
        ylabel('WPE'); 
        text(3.325,1.2,'C)','FontSize', 12, 'FontWeight', 'Bold');

subplot(4,1,4)
w=5;
        h3=scatter(WPE_data(WPE_data(1:end-4,2)==t,w),RMSE(1:end-1),30,WPE_data(WPE_data(1:end-4,2)==t,1),'o','filled');
        colormap;
 
            set(gca, ...  'yticklabel',{'0','0.1','0.2','0.3','0.4'} , ...
            'YLim', [-0.00003, 0.00071]    , ... 
            'YTick', 0:0.0002:0.0006      , ...
            'yticklabel',{'0','2','4','6'} , ...
            'XLim', [0.075, 0.65]    , ...
            'XTick', 0.1:0.1:0.6      , ...
            'xticklabel',{'0.1','0.2','0.3','0.4','0.5','0.6'} , ...
            'Linewidth', 1, ...
            'Box','on'     ) 
        
        ArStart=[0.09 0.00119; 0.45 0.0014; 0.58 0.00123];
        
        for ii=1:length(r_plot)
            iiy=find(round(rs,10)==r_plot(ii),1); 
            iix=find(round(WPE_data(:,1),10)==r_plot(ii) & WPE_data(:,2)==t);
            arrow([ArStart(ii,1),ArStart(ii,2)], [WPE_data(iix,w),RMSE(iiy)+0.00002], 'Color', arcol(ii,:),'BaseAngle',60)
        end
               
        
        ylabel('nRMSE')  
                xlabel('WPE')
        text(-0.01,0.0008,'D)','FontSize', 12, 'FontWeight', 'Bold');
        text(0.01,0.0007,'x10^-^4','FontSize', 12);
        
saveas(f4,'Fig_4_Logistic_map_Bifurc_WPE_RMSE','pdf')        



%% plot nRMSE as a function of PE for different levels of noise in separate subplots  (Fig 5 in MS)  

RMSE_Zpn=output.rmse(strcmp(output.noise_level,'0_zero') & strcmp(output.noise_type,'process_noise'));
RMSE_Lpn=output.rmse(strcmp(output.noise_level,'1_low') & strcmp(output.noise_type,'process_noise'));
RMSE_Mpn=output.rmse(strcmp(output.noise_level,'2_medium') & strcmp(output.noise_type,'process_noise'));
RMSE_Hpn=output.rmse(strcmp(output.noise_level,'3_high') & strcmp(output.noise_type,'process_noise'));
RMSE_Zon=output.rmse(strcmp(output.noise_level,'0_zero') & strcmp(output.noise_type,'obs_noise'));
RMSE_Lon=output.rmse(strcmp(output.noise_level,'1_low') & strcmp(output.noise_type,'obs_noise'));
RMSE_Mon=output.rmse(strcmp(output.noise_level,'2_medium') & strcmp(output.noise_type,'obs_noise'));
RMSE_Hon=output.rmse(strcmp(output.noise_level,'3_high') & strcmp(output.noise_type,'obs_noise'));

%Redefine these variables otherwise code returns an error as it thinks they might be functions... IDK
WPE_Zpn=WPE_Zpn;
WPE_Lpn=WPE_Lpn;
WPE_Mpn=WPE_Mpn;
WPE_Hpn=WPE_Hpn;
WPE_Zon=WPE_Zon;
WPE_Lon=WPE_Lon;
WPE_Mon=WPE_Mon;
WPE_Hon=WPE_Hon;

f5=figure(5);clf
hold all
set(f5, 'Position', [0, 0, 800, 400]);
set(f5, 'DefaultAxesFontSize', 14);
set(f5, 'DefaultLineLineWidth',1.5);
set(f5, 'color','white')    

w=5; %WPE word length to plot
t=1; %WPE tau to plot   

%Process noise
sp1=subplot(3,2,1);
scatter(WPE_Zpn(WPE_Zpn(1:end-16,2)==t,w),RMSE_Zpn(1:end-4),30,[0.7,0.7,0.7],'o'); hold on; 
hh4=scatter(WPE_Hpn(WPE_Hpn(1:end-16,2)==t,w),RMSE_Hpn(1:end-4),30,WPE_Hpn(WPE_Hpn(1:end-16,2)==t,1),'o','filled'); hold on; 
    text(0.275,0.021,'A) Process noise','FontSize', 14, 'FontWeight', 'Bold');        
    text(0.3,0.0146,'High, SD = 0.01','FontSize', 14);        

  set(gca, ...  
            'YLim', [-0.001, 0.018]    , ... 
            'YTick', [0, 0.008, 0.016]      , ...
            'yticklabel',{'0', '0.008','0.016'} , ...
            'XLim', [0.275, 0.725]    , ...
            'XTick', 0.3:0.1:0.7      , ...
            'xticklabel',{} , ...
            'Linewidth', 1, ...
            'Box','on'     ) 
        
sp3=subplot(3,2,3);
scatter(WPE_Zpn(WPE_Zpn(1:end-16,2)==t,w),RMSE_Zpn(1:end-4),30,[0.7,0.7,0.7],'o'); hold on; 
hm4=scatter(WPE_Mpn(WPE_Mpn(1:end-16,2)==t,w),RMSE_Mpn(1:end-4),30,WPE_Mpn(WPE_Mpn(1:end-16,2)==t,1),'o','filled'); hold on; 
    ylabel('nRMSE')
    text(0.3,0.0046,'Med, SD = 0.001','FontSize', 14);   
    
  set(gca, ...  
            'YLim', [-0.0004, 0.0055]    , ... 
            'YTick', [0, 0.002, 0.004]      , ...
            'yticklabel',{'0', '0.002','0.004'} , ...
            'XLim', [0.275, 0.725]    , ...
            'XTick', 0.3:0.1:0.7      , ...
            'xticklabel',{} , ...
            'Linewidth', 1, ...
            'Box','on'     )
       
 pos=get(sp3,'Position');   set(sp3, 'Position',[pos(1) pos(2)+0.08 pos(3) pos(4)])
                
sp5=subplot(3,2,5);
scatter(WPE_Zpn(WPE_Zpn(1:end-16,2)==t,w),RMSE_Zpn(1:end-4),30,[0.7,0.7,0.7],'o'); hold on; 
hl4=scatter(WPE_Lpn(WPE_Lpn(1:end-16,2)==t,w),RMSE_Lpn(1:end-4),30,WPE_Lpn(WPE_Lpn(1:end-16,2)==t,1),'o','filled'); hold on; 
    xlabel('WPE')
    text(0.3,0.0018,'Low, SD = 0.0001','FontSize', 14);   
  set(gca, ...  
            'YLim', [-0.00015, 0.0022]    , ... 
            'YTick', [0, 0.001, 0.002]      , ...
            'yticklabel',{'0', '0.001','0.002'} , ...
            'XLim', [0.275, 0.725]    , ...
            'XTick', 0.3:0.1:0.7      , .....
            'xticklabel',{'0.3','0.4','0.5','0.6','0.7'} , ...
            'Linewidth', 1, ...
            'Box','on'     )

 pos=get(sp5,'Position');   set(sp5, 'Position',[pos(1) pos(2)+0.16 pos(3) pos(4)])


 
%Observational noise
sp2=subplot(3,2,2);
scatter(WPE_Zon(WPE_Zon(1:end-16,2)==t,w),RMSE_Zon(1:end-4),30,[0.7,0.7,0.7],'o'); hold on; 
hh4=scatter(WPE_Hon(WPE_Hon(1:end-16,2)==t,w),RMSE_Hon(1:end-4),30,WPE_Hon(WPE_Hon(1:end-16,2)==t,1),'o','filled'); hold on; 
    text(0.275,0.021,'B) Observational noise','FontSize', 14, 'FontWeight', 'Bold');        

  set(gca, ...  
           'YLim', [-0.001, 0.018]    , ... 
            'YTick', [0, 0.008, 0.016]      , ...
            'yticklabel',{} , ...
            'XLim', [0.275, 0.725]    , ...
            'XTick', 0.3:0.1:0.7      , ...
            'xticklabel',{} , ...
            'Linewidth', 1, ...
            'Box','on'     ) 

pos=get(sp2,'Position');   set(sp2, 'Position',[pos(1)-0.1 pos(2) pos(3) pos(4)])        
        
sp4=subplot(3,2,4) ;
scatter(WPE_Zon(WPE_Zon(1:end-16,2)==t,w),RMSE_Zon(1:end-4),30,[0.7,0.7,0.7],'o'); hold on; 
hm4=scatter(WPE_Mon(WPE_Mon(1:end-16,2)==t,w),RMSE_Mon(1:end-4),30,WPE_Mon(WPE_Mon(1:end-16,2)==t,1),'o','filled'); hold on; 
    
  set(gca, ...  
            'YLim', [-0.0004, 0.0055]    , ... 
            'YTick', [0, 0.002, 0.004]      , ...
            'yticklabel',{} , ...
            'XLim', [0.275, 0.725]    , ...
            'XTick', 0.3:0.1:0.7      , ...
            'xticklabel',{} , ...
            'Linewidth', 1, ...
            'Box','on'     )
       
 pos=get(sp4,'Position');   set(sp4, 'Position',[pos(1)-0.1 pos(2)+0.08 pos(3) pos(4)])
 
sp6=subplot(3,2,6);
scatter(WPE_Zon(WPE_Zon(1:end-16,2)==t,w),RMSE_Zon(1:end-4),30,[0.7,0.7,0.7],'o'); hold on; 
hl4=scatter(WPE_Lon(WPE_Lon(1:end-16,2)==t,w),RMSE_Lon(1:end-4),30,WPE_Lon(WPE_Lon(1:end-16,2)==t,1),'o','filled'); hold on; 
    xlabel('WPE')
  set(gca, ...  
            'YLim', [-0.00015, 0.0022]    , ... 
            'YTick', [0, 0.001, 0.002]      , ...
            'yticklabel',{} , ...
            'XLim', [0.275, 0.725]    , ...
            'XTick', 0.3:0.1:0.7      , ...
            'xticklabel',{'0.3','0.4','0.5','0.6','0.7'} , ...
            'Linewidth', 1, ...
            'Box','on'     )
        
 pos=get(sp6,'Position');   set(sp6, 'Position',[pos(1)-0.1 pos(2)+0.16 pos(3) pos(4)])
 
 orient(f5,'landscape')
print(f5,'Fig_5_Logistic_map_WPE_RMSE_noise','-dpdf')

