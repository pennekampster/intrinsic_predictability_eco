%Calculate the WPE of the Ricker Model with lognormal, multiplicative noise across various levels of Tau and
%Permutation length

clearvars
close all

%% Simulate the Ricker model   N(t+1)=N(t)exp(r*(1-n(t)/K))
% with lognormal process and observational noise, calculate WPE, save data file for analysis of EDM in R

M = 10000; % number of iterations of ricker model
N = 3000; % length of end of time series to save and analyzes

R = 500; % number of "r" values to simulate
rmin = 2.4; % starting value of "r"
rmax = 3.4; % final value of "r"... anything higher diverges.
    rs = linspace(rmin,rmax,R); % vector of "r" values
K=1000; %Carrying capacity

mwl=5; % Maximum word length of permutations to consider
mtau=4; % Maximum time delay between the 'letters' of each 'word' to consider

scale = 10000; % determines the level of rounding for plotting the bifurcation diagram (not for analysis)
maxpoints = 200; % determines maximum values to plot

pn = [0;0.001;0.01;0.1]; % standard deviation of the process noise; levels of process noise to loop through
RM_WPE_Zpn=[]; RM_WPE_Lpn=[]; RM_WPE_Mpn=[]; RM_WPE_Hpn=[];% allocate matrices to store WPE output for zero, low, medium and high process noise
RM_Zpn=zeros(N,R); RM_Lpn=zeros(N,R); RM_Mpn=zeros(N,R); RM_Hpn=zeros(N,R); % allocate matrices to store LM time series output for zero, low, medium and high process noise
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
    
    for i = 2:M % iterate the Ricker model with process noise
                
        %Add process noise to x at every time step, analogous to migration
        %becasuse the growth rate is not changing, just the population
        %size.   
        x(i) = x(i-1)*exp(r*(1-(x(i-1)/K)))*exp(randn(1,1).*pn(k)-((pn(k)^2)/2)); %Ricker model with lognormal process noise. See Hilborn & Mangel 1997 The Ecological Detective, pg 146 to justify this formulation of the lognormal process error

         x_expct=x(i-1)*exp(r*(1-(x(i-1)/K))); %the expected x with no noise
         xdiff(i)=x(i)-x_expct; %the difference between x and the expected x
        
    end
    
    % Calculate the Weighted Permutation Entropy for Ricker with process noise
        % Loop through the values of tau
        x(1:M-N,:)=[]; %Remove initial transient dynamics
        
        xr=round(x, 10); %Round to nearest 10^-10 to remove minor rounding fluctuations
        for t=1:mtau
            tau=t;
        [ H ] = WPE( xr, mwl, tau );
        
            if k==1; RM_WPE_Zpn=[RM_WPE_Zpn; r, tau, H, pn(k)];
            elseif k==2; RM_WPE_Lpn=[RM_WPE_Lpn; r, tau, H, pn(k)];
            elseif k==3; RM_WPE_Mpn=[RM_WPE_Mpn; r, tau, H, pn(k)];
            elseif k==4; RM_WPE_Hpn=[RM_WPE_Hpn; r, tau, H, pn(k)];
            end
        end
    
    % only save those unique, semi-stable values from zero process
    % noise to plot in bifurcation plot.
        if k==1; out{j} = unique(round(scale*x(end-maxpoints:end))); end % only save bifurcation data for zero added noise otherwise the plot is messy and unrecognizable
    % save LM time series data for export to R
        if k==1; RM_Zpn(:,j)=xr; xdiff_Zpn(:,j)=xdiff;
        elseif k==2; RM_Lpn(:,j)=xr; xdiff_Lpn(:,j)=xdiff;
        elseif k==3; RM_Mpn(:,j)=xr; xdiff_Mpn(:,j)=xdiff;
        elseif k==4; RM_Hpn(:,j)=xr; xdiff_Hpn(:,j)=xdiff;
        end
        
        
end %r values loop

end %process noise loop

% Adding observational noise
    on = [0; 0.001; 0.01; 0.1]; % standard deviation of the observational noise levels to loop through

    % Multiply observations by lognormally distributed error. See Hilborn
    % and Mangel 1997 The Ecological Detective for an explanation of this
    % formulation of lognormal error
    RM_Zon=exp(randn(N,R).*on(1)-((on(1)^2)/2)).*RM_Zpn;
    RM_Lon=exp(randn(N,R).*on(2)-((on(2)^2)/2)).*RM_Zpn;
    RM_Mon=exp(randn(N,R).*on(3)-((on(3)^2)/2)).*RM_Zpn;
    RM_Hon=exp(randn(N,R).*on(4)-((on(4)^2)/2)).*RM_Zpn;
    
    % Calculate the Weighted Permutation Entropy for observational noise
      RM_WPE_Zon=[]; RM_WPE_Lon=[]; RM_WPE_Mon=[]; RM_WPE_Hon=[];  % allocate matrices to store WPE output for zero, low, medium, high and superhigh observational noise

            for j = 1:length(rs) % Loop through the different r values
                r=rs(j); 
                
            for t=1:mtau % Loop through the values of tau for WPE calculation
                tau=t;    
                   [ H ] = WPE(RM_Zon(:,j), mwl, tau);
                   RM_WPE_Zon=[RM_WPE_Zon; r, tau, H, on(1)];

                   [ H ] = WPE(RM_Lon(:,j), mwl, tau);
                   RM_WPE_Lon=[RM_WPE_Lon; r, tau, H, on(2)];

                   [ H ] = WPE(RM_Mon(:,j), mwl, tau);
                   RM_WPE_Mon=[RM_WPE_Mon; r, tau, H, on(3)];
%
                   [ H ] = WPE(RM_Hon(:,j), mwl, tau);
                   RM_WPE_Hon=[RM_WPE_Hon; r, tau, H, on(4)];

            end %tau loop
            end %r loop    
        
save 'RickerTimeSeriesData.mat' RM_Zpn RM_Lpn RM_Mpn RM_Hpn RM_Zon RM_Lon RM_Mon RM_Hon rs; %Save time series so to access from R and calcuate the S-map projection and the absolute scales error
save 'RickerWPEData.mat' RM_WPE_Zpn RM_WPE_Lpn RM_WPE_Mpn RM_WPE_Hpn RM_WPE_Zon RM_WPE_Lon RM_WPE_Mon RM_WPE_Hon; %Save time series so to access from R and calcuate the S-map projection and the absolute scales error



%% Figure WPE vs r

% Rearrange cell array into a large n-by-2 vector for plotting
bfrctn = [];
for k = 1:length(rs)
    n = length(out{k});
    bfrctn = [bfrctn;  rs(k)*ones(n,1),out{k}];
end

% Plot the data
WPE_data=RM_WPE_Zpn;% plot the zero process noise results
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
            'XLim', [rmin-0.1, rmax+0.1]         , ...
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

%% Plot of the mean absolute scaled error from S-Map projections
% First do this: load 'RickerTimeSeriesData.mat' into the R script
% 'EDM_Ricker_noise_for_Matlab_v2.R, which uses the function 'ASE.R' to calculate the
% root mean squared error of the Ricker timeseries at each of 100 
% projected time steps using the previous 1900 timestep to formulate the 
% attractor of the S-Map. ASE is just zero when at steady state 
% equilibrium (r<3) or when the population crashes (r>=4).

load 'RickerTimeSeriesData.mat' 
load 'RickerWPEData.mat'

mwl=5; 
mtau=4; 
N=3000;

load Ricker_Forecast_Error.mat; %load nRMSE results of simulations of Ricker model with zero process noies
RMSE=output.rmse(strcmp(output.noise_level,'0_zero') & strcmp(output.noise_type,'process_noise'));
WPE_data=RM_WPE_Zpn;% plot the zero process noise results

f2=figure(2);clf
hold all
set(f2, 'Position', [100, 100, 1000, 1000]);
set(f2, 'DefaultAxesFontSize', 12);
set(f2, 'DefaultLineLineWidth',1.2);
set(f2, 'color','white')
Alphabet = 'abcdefghijkl';

for t=1:mtau
    for w=3:mwl
        p=((t-1)*3)+w-2;
        sp=subplot(4,3,p);
        
        h2=scatter(WPE_data(WPE_data(:,2)==t,w),RMSE,30,WPE_data(WPE_data(:,2)==t,1));
        colormap jet;
        cm=colormap;
        cm=cm(1:64,:);
        colormap(cm);
        hold on
        
        title(strcat(Alphabet(p),') \tau=',num2str(t),', wl=', num2str(w)));
        
        set(gca, ... 
            'XTick', 0:0.2:1           , ...
            'XLim', [0, 1]         , ... %'YTick', 0:0.004:0.02      , ... %'YLim', [0, 0.02]    , ...
            'xticklabel',[] , ...
            'yticklabel',[] , ...
            'Box', 'on'             , ...
            'Linewidth', 1          , ...
            'TickDir','out')

        if p==1 || p==4 || p==7 || p==10, set(gca, 'yticklabel',{'0.0','0.004','0.008','0.012','0.016','0.02'}),end
        if p==10 || p==11 || p==12 , set(gca, 'xticklabel',{'0','0.2','0.4','0.6','0.8','1'}),end

        %Modify the positions/size of each subplot so that there is less white space
        pos = get(sp, 'position'); %Get positions of each subplot
        if p==2 || p==5 || p==8 || p==11,  pos(1)=pos(1)+0.03; end
        if p==1 || p==4 || p==7 || p==10,   pos(1)=pos(1)+0.06; end
        if p==4 || p==5 || p==6,   pos(2)=pos(2)+0.03;  end
        if p==7 || p==8 || p==9,    pos(2)=pos(2)+0.06; end
        if p==10 || p==11 || p==12,   pos(2)=pos(2)+0.09;  end
        set(sp, 'position', pos);
        
        if p==11, xlabel('WPE'); end
        if p==7, yL=ylabel('RMSE'); end
        
 

        if p==3, h2cb=colorbar('YTick',1.5:0.5:4);
        set(get(h2cb,'Title'),'String','r'); end


    end
end


%% Main RM Summary figure

f3=figure(3);clf
hold all
set(f3, 'Position', [100, 100, 400, 800]);
set(f3, 'DefaultAxesFontSize', 14);
set(f3, 'DefaultLineLineWidth',1.5);
set(f3, 'color','white')

t=1; %The chosen tau to plot

subplot(4,1,1) %Demonstrate the conceptual leap from time series to bifurcation diagram. Plot 3 different time series of 3 different growth rates in the corresponding colors.

% Choose growth rate values to highlight in figures
r_plot=[2.4200400802 2.8529058116 3.1434869739];

col=colormap;

xvals=linspace(2970,3000,30)';

arcol=zeros(length(r_plot),3);

for l=1:length(r_plot)
    ll=find(round(rs,10)==r_plot(l),1); 
    yvals=RM_Zpn(2971:3000,ll);
    curcol=round(((length(col)-1)/length(rs))*ll)+1; %point color on plot
    arcol(l,:)=col(curcol,:);
    h1=plot(xvals,yvals, 'o-', 'Color',col(curcol,:), 'Linewidth', 1.5, 'MarkerFaceColor',col(curcol,:));
        hold on   
end

 set(gca, ... 
            'YLim', [-100, 3000]    , ...
            'YTick', 0:1000:3000      , ...
            'yticklabel',{'0','1000','2000','3000'} , ...
            'XLim', [N-25, N]         , ...
            'Linewidth', 1, ...
            'Box','on')
           
        xlabel('time step')
        ylabel('Abundance'); 
        text(2975,3400,'A)','FontSize', 12, 'FontWeight', 'Bold');


subplot(4,1,2)
        h3=scatter(bfrctn(1:end-1,1),bfrctn(1:end-1,2)/scale,2,bfrctn(1:end-1,1),'o','filled');
        hold on
       
        set(gca, 'XLim', [rmin-0.01, rmax+0.01]         , ...
            'XTick', 2.4:0.2:3.4      , ...
            'xticklabel',{'2.4','2.6','2.8','3.0','3.2','3.4'} , ...
            'Linewidth', 1, ...
            'YLim', [-100, 3500]    , ...
            'YTick', 0:1000:3000      , ...
            'yticklabel',{'0','1000','2000','3000'} , ...
            'Box','on')
        
        
         for ii=1:length(r_plot)
            plot([r_plot(ii),r_plot(ii)],[3600,-100], '.-', 'Color',arcol(ii,:), 'Linewidth', 2);
        hold on   
         end 
         
        xlabel('Growth rate, r')
        ylabel('Abundance'); 
        text(2.39,3800,'B)','FontSize', 12, 'FontWeight', 'Bold');

subplot(4,1,3) 
WPE_data=RM_WPE_Zpn;% plot the zero process noise results

plot(WPE_data(WPE_data(1:end-4,2)==t,1),WPE_data(WPE_data(1:end-4,2)==t,3),'-', 'Color',[0.7,0.7,0.7],'Linewidth', 1.5); hold on
plot(WPE_data(WPE_data(1:end-4,2)==t,1),WPE_data(WPE_data(1:end-4,2)==t,4),'-', 'Color',[0.4,0.4,0.4],'Linewidth', 1.5); hold on
plot(WPE_data(WPE_data(1:end-4,2)==t,1),WPE_data(WPE_data(1:end-4,2)==t,5),'k-', 'Linewidth', 1.5); hold on

      set(gca, ... 
            'XLim', [rmin-0.01, rmax+0.01]         , ...
            'XTick', 2.4:0.2:3.4      , ...
            'xticklabel',{'2.4','2.6','2.8','3.0','3.2','3.4'} , ...
            'Linewidth', 1, ...
            'Box','on')
     
        for ii=1:length(r_plot)
            iiy=find(round(rs,10)==r_plot(ii),1); 
            iix=find(round(WPE_data(:,1),10)==r_plot(ii) & WPE_data(:,2)==t);
            scatter(r_plot(ii),WPE_data(iix,5), 60, arcol(ii,:),'o','filled')
        end
        
        xlabel('Growth rate, r')
        ylabel('WPE'); 
        text(2.39,1.1,'C)','FontSize', 12, 'FontWeight', 'Bold');

subplot(4,1,4)
w=5;
        h3=scatter(WPE_data(WPE_data(1:end-4,2)==t,w),RMSE(1:end-1),30,WPE_data(WPE_data(1:end-4,2)==t,1),'o','filled');
        colormap;
 
            set(gca, ...  
            'YLim', [-10, 150]    , ...
            'YTick', 0:50:150     , ...
            'yticklabel',{'0','50','100', '150'} , ...   
            'XLim', [0.1, 0.6]    , ...
            'XTick', 0.1:0.1:0.6      , ...
            'xticklabel',{'0.1','0.2','0.3','0.4','0.5','0.6'} , ...
            'Linewidth', 1, ...
            'Box','on'     ) 
        
        ArStart=[0.125 250; 0.33 300; 0.475 260];
        
        for ii=1:length(r_plot)
            iiy=find(round(rs,10)==r_plot(ii),1); 
            iix=find(round(WPE_data(:,1),10)==r_plot(ii) & WPE_data(:,2)==t);
            arrow([ArStart(ii,1),ArStart(ii,2)], [WPE_data(iix,w),RMSE(iiy)+0.00002], 'Color', arcol(ii,:),'BaseAngle',60)
        end
               
        
        ylabel('nRMSE')  
        xlabel('WPE')
        text(0.1,165,'D)','FontSize', 12, 'FontWeight', 'Bold');
        
      
print(f3,'Fig_Ricker_Bifurc_WPE_RMSE','-dpdf')


%% plot nRMSE as a function of PE for different levels of process noise, rainbow colors connected by lines       

RMSE_Zpn=output.rmse(strcmp(output.noise_level,'0_zero') & strcmp(output.noise_type,'process_noise'));
RMSE_Lpn=output.rmse(strcmp(output.noise_level,'1_low') & strcmp(output.noise_type,'process_noise'));
RMSE_Mpn=output.rmse(strcmp(output.noise_level,'2_medium') & strcmp(output.noise_type,'process_noise'));
RMSE_Hpn=output.rmse(strcmp(output.noise_level,'3_high') & strcmp(output.noise_type,'process_noise'));
RMSE_Zon=output.rmse(strcmp(output.noise_level,'0_zero') & strcmp(output.noise_type,'obs_noise'));
RMSE_Lon=output.rmse(strcmp(output.noise_level,'1_low') & strcmp(output.noise_type,'obs_noise'));
RMSE_Mon=output.rmse(strcmp(output.noise_level,'2_medium') & strcmp(output.noise_type,'obs_noise'));
RMSE_Hon=output.rmse(strcmp(output.noise_level,'3_high') & strcmp(output.noise_type,'obs_noise'));

f4=figure(4);clf
hold all
set(f4, 'Position', [100, 100, 600, 900]);
set(f4, 'DefaultAxesFontSize', 14);
set(f4, 'DefaultLineLineWidth',1);
set(f4, 'color','white')    

w=5;
t=1;    
col=colormap;

subplot(2,1,1) %Process noise

XWPE = horzcat(RM_WPE_Zpn(RM_WPE_Zpn(:,2)==t,w), RM_WPE_Lpn(RM_WPE_Lpn(:,2)==t,w), RM_WPE_Mpn(RM_WPE_Mpn(:,2)==t,w), RM_WPE_Hpn(RM_WPE_Hpn(:,2)==t,w));
YnRMSE = horzcat(RMSE_Zpn, RMSE_Lpn, RMSE_Mpn, RMSE_Hpn);
ZGR = horzcat(rs',rs',rs',rs'); %Growth rate on the Z axis

for ii=1:length(XWPE)
        curcol=round(((length(col)-1)/length(XWPE))*ii)+1; %current color
        plot3(XWPE(ii,:),ZGR(ii,:),YnRMSE(ii,:),'Color',col(curcol,:)); hold on;
end
view(3)
plot3(XWPE(:,1),ZGR(:,1),YnRMSE(:,1),'-','Color',[0 0 0]); hold on; %add black no-noise baseline


        zlabel('Residual mean squared error, RMSE (x10^-^3)','FontSize', 14)
        xlabel('WPE','FontSize', 14)
        ylabel('Growth rate, r','FontSize', 14)
        title('A) Process noise','FontSize', 14, 'FontWeight', 'Bold');

        xh = get(gca,'XLabel'); % Handle of the x label
        set(xh, 'Units', 'Normalized')
        pos = get(xh, 'Position');
        set(xh, 'Position',pos.*[1.1,.5,0],'Rotation',12)
        yh = get(gca,'YLabel'); % Handle of the y label
        set(yh, 'Units', 'Normalized')
        pos = get(yh, 'Position');
        set(yh, 'Position',pos.*[1,3,1],'Rotation',-24)
        zh = get(gca,'ZLabel'); % Handle of the z label
        set(zh, 'Units', 'Normalized')
        pos = get(zh, 'Position');
        set(zh, 'Position',pos.*[1,1,0],'Rotation',90)
        
        
            set(gca, ... 
            'Ydir','reverse'    , ...
            'YLim', [1.8, 3.4]    , ...
            'YTick', 1.8:0.4:3.4      , ...
            'yticklabel',{'1.8','2.2','2.6','3.0','3.4'} , ...
            'Linewidth', 1 , ...
            'Box', 'off')   
        
subplot(2,1,2) %Observational noise

XWPE = horzcat(RM_WPE_Zon(RM_WPE_Zon(:,2)==t,w), RM_WPE_Lon(RM_WPE_Lon(:,2)==t,w), RM_WPE_Mon(RM_WPE_Mon(:,2)==t,w), RM_WPE_Hon(RM_WPE_Hon(:,2)==t,w));
YnRMSE = horzcat(RMSE_Zon, RMSE_Lon, RMSE_Mon, RMSE_Hon);
ZGR = horzcat(rs',rs',rs',rs'); %Growth rate on the Z axis

for ii=1:length(XWPE)
        curcol=round(((length(col)-1)/length(XWPE))*ii)+1; %current color
        plot3(XWPE(ii,:),ZGR(ii,:),YnRMSE(ii,:),'Color',col(curcol,:)); hold on;
end
view(3)
plot3(XWPE(:,1),ZGR(:,1),YnRMSE(:,1),'-','Color',[0 0 0]); hold on; %add black no-noise baseline


         zlabel('Residual mean squared error, RMSE (x10^-^3)','FontSize', 14)
        xlabel('WPE','FontSize', 14)
        ylabel('Growth rate, r','FontSize', 14)
        title('B) Observational noise','FontSize', 14, 'FontWeight', 'Bold');
        
 
        xh = get(gca,'XLabel'); % Handle of the x label
        set(xh, 'Units', 'Normalized')
        pos = get(xh, 'Position');
        set(xh, 'Position',pos.*[1.1,0,0],'Rotation',12)
        yh = get(gca,'YLabel'); % Handle of the y label
        set(yh, 'Units', 'Normalized')
        pos = get(yh, 'Position');
        set(yh, 'Position',pos.*[1,2,1],'Rotation',-24)
        zh = get(gca,'ZLabel'); % Handle of the z label
        set(zh, 'Units', 'Normalized')
        pos = get(zh, 'Position');
        set(zh, 'Position',pos.*[1,1,0],'Rotation',90)
        
        
            set(gca, ... 
            'Ydir','reverse'    , ...
            'YLim', [1.8, 3.4]    , ...
            'YTick', 1.8:0.4:3.4      , ...
            'yticklabel',{'1.8','2.2','2.6','3.0','3.4'} , ...
            'Linewidth', 1 , ...
            'Box', 'off')   
        
        
        
       

%% plot nRMSE as a function of PE for different levels of noise in separate subplots     
f5=figure(5);clf
hold all
set(f5, 'Position', [0, 0, 800, 400]);
set(f5, 'DefaultAxesFontSize', 14);
set(f5, 'DefaultLineLineWidth',1.5);
set(f5, 'color','white')    

w=5; %WPE word length to plot
t=1; %WPE tau to plot   

%Weirdly need to redefine these as variable? Otherwise code returns an
%error as it thinks they might be functions... IDK
RM_WPE_Zpn=RM_WPE_Zpn;
RM_WPE_Lpn=RM_WPE_Lpn;
RM_WPE_Mpn=RM_WPE_Mpn;
RM_WPE_Hpn=RM_WPE_Hpn;
RM_WPE_Zon=RM_WPE_Zon;
RM_WPE_Lon=RM_WPE_Lon;
RM_WPE_Mon=RM_WPE_Mon;
RM_WPE_Hon=RM_WPE_Hon;



%Process noise
sp1=subplot(3,2,1);
scatter(RM_WPE_Zpn(RM_WPE_Zpn(1:end-16,2)==t,w),RMSE_Zpn(1:end-4),30,[0.7,0.7,0.7],'o'); hold on; 
hh4=scatter(RM_WPE_Hpn(RM_WPE_Hpn(1:end-16,2)==t,w),RMSE_Hpn(1:end-4),30,RM_WPE_Hpn(RM_WPE_Hpn(1:end-16,2)==t,1),'o','filled'); hold on; 
    text(0.1,280,'A) Process noise','FontSize', 14, 'FontWeight', 'Bold');        
    text(0.12,220, horzcat('High, SD = ', num2str(pn(4))), 'FontSize', 14);        

    set(gca, ...  
            'YLim', [-10, 250]    , ...
            'YTick', 0:100:200     , ...
            'yticklabel',{'0','100','200'} , ...   
            'XLim', [0.1, 0.7]    , ...
            'XTick', 0.1:0.1:0.7      , ...
            'xticklabel',{} , ...
            'Linewidth', 1, ...
            'Box','on'     ) 


        
sp3=subplot(3,2,3);
scatter(RM_WPE_Zpn(RM_WPE_Zpn(1:end-16,2)==t,w),RMSE_Zpn(1:end-4),30,[0.7,0.7,0.7],'o'); hold on; 
hm4=scatter(RM_WPE_Mpn(RM_WPE_Mpn(1:end-16,2)==t,w),RMSE_Mpn(1:end-4),30,RM_WPE_Mpn(RM_WPE_Mpn(1:end-16,2)==t,1),'o','filled'); hold on; 
    ylabel('nRMSE')
    text(0.12,220, horzcat('Med, SD = ', num2str(pn(3))),'FontSize', 14);   
    
    set(gca, ...  
            'YLim', [-10, 250]    , ...
            'YTick', 0:100:200     , ...
            'yticklabel',{'0','100','200'} , ...
            'XLim', [0.1, 0.7]    , ...
            'XTick', 0.1:0.1:0.7      , ...
            'xticklabel',{} , ...
            'Linewidth', 1, ...
            'Box','on'     ) 

       
 pos=get(sp3,'Position');   set(sp3, 'Position',[pos(1) pos(2)+0.08 pos(3) pos(4)])
                
sp5=subplot(3,2,5);
scatter(RM_WPE_Zpn(RM_WPE_Zpn(1:end-16,2)==t,w),RMSE_Zpn(1:end-4),30,[0.7,0.7,0.7],'o'); hold on; 
hl4=scatter(RM_WPE_Lpn(RM_WPE_Lpn(1:end-16,2)==t,w),RMSE_Lpn(1:end-4),30,RM_WPE_Lpn(RM_WPE_Lpn(1:end-16,2)==t,1),'o','filled'); hold on; 
    xlabel('WPE')
    text(0.12,220, horzcat('Low, SD = ', num2str(pn(2))),'FontSize', 14);   
    set(gca, ...  
            'YLim', [-10, 250]    , ...
            'YTick', 0:100:200     , ...
            'yticklabel',{'0','100','200'} , ...
            'XLim', [0.1, 0.7]    , ...
            'XTick', 0.1:0.1:0.7      , ...
            'xticklabel',{'0.1','0.2','0.3','0.4','0.5','0.6',''} , ...
            'Linewidth', 1, ...
            'Box','on'     ) 


 pos=get(sp5,'Position');   set(sp5, 'Position',[pos(1) pos(2)+0.16 pos(3) pos(4)])


 
%Observational noise
sp2=subplot(3,2,2);
scatter(RM_WPE_Zon(RM_WPE_Zon(1:end-16,2)==t,w),RMSE_Zon(1:end-4),30,[0.7,0.7,0.7],'o'); hold on; 
hh4=scatter(RM_WPE_Hon(RM_WPE_Hon(1:end-16,2)==t,w),RMSE_Hon(1:end-4),30,RM_WPE_Hon(RM_WPE_Hon(1:end-16,2)==t,1),'o','filled'); hold on; 
    text(0.1,280,'B) Observational noise','FontSize', 14, 'FontWeight', 'Bold');        

    set(gca, ...  
            'YLim', [-10, 250]    , ...
            'YTick', 0:100:200     , ...
            'yticklabel',{} , ...   
            'XLim', [0.1, 0.7]    , ...
            'XTick', 0.1:0.1:0.7      , ...
            'xticklabel',{} , ...
            'Linewidth', 1, ...
            'Box','on'     ) 

pos=get(sp2,'Position');   set(sp2, 'Position',[pos(1)-0.1 pos(2) pos(3) pos(4)])        
        
sp4=subplot(3,2,4) ;
scatter(RM_WPE_Zon(RM_WPE_Zon(1:end-16,2)==t,w),RMSE_Zon(1:end-4),30,[0.7,0.7,0.7],'o'); hold on; 
hm4=scatter(RM_WPE_Mon(RM_WPE_Mon(1:end-16,2)==t,w),RMSE_Mon(1:end-4),30,RM_WPE_Mon(RM_WPE_Mon(1:end-16,2)==t,1),'o','filled'); hold on; 
    
    set(gca, ...  
            'YLim', [-10, 250]    , ...
            'YTick', 0:100:200     , ...
             'yticklabel',{} , ...   
            'XLim', [0.1, 0.7]    , ...
            'XTick', 0.1:0.1:0.7      , ...
            'xticklabel',{} , ...
            'Linewidth', 1, ...
            'Box','on'     ) 
        
 pos=get(sp4,'Position');   set(sp4, 'Position',[pos(1)-0.1 pos(2)+0.08 pos(3) pos(4)])
 
sp6=subplot(3,2,6);
scatter(RM_WPE_Zon(RM_WPE_Zon(1:end-16,2)==t,w),RMSE_Zon(1:end-4),30,[0.7,0.7,0.7],'o'); hold on; 
hl4=scatter(RM_WPE_Lon(RM_WPE_Lon(1:end-16,2)==t,w),RMSE_Lon(1:end-4),30,RM_WPE_Lon(RM_WPE_Lon(1:end-16,2)==t,1),'o','filled'); hold on; 
    xlabel('WPE')
    set(gca, ...  
 'YLim', [-10, 250]    , ...
            'YTick', 0:100:200     , ...
            'yticklabel',{} , ...
             'XLim', [0.1, 0.7]    , ...
            'XTick', 0.1:0.1:0.7      , ...
            'xticklabel',{'0.1','0.2','0.3','0.4','0.5','0.6',''} , ...
            'Linewidth', 1, ...
            'Box','on'     ) 
        
 pos=get(sp6,'Position');   set(sp6, 'Position',[pos(1)-0.1 pos(2)+0.16 pos(3) pos(4)])
 
 orient(f5,'landscape')
print(f5,'Fig_Ricker_WPE_RMSE_noise_sep_subplots','-dpdf','-bestfit')


%% plot RMSE as a function of PE for different levels of noise all together in one plot
%Weirdly need to redefine these as variable? Otherwise code returns an
%error as it thinks they might be functions... IDK
RM_WPE_Zpn=RM_WPE_Zpn;
RM_WPE_Lpn=RM_WPE_Lpn;
RM_WPE_Mpn=RM_WPE_Mpn;
RM_WPE_Hpn=RM_WPE_Hpn;
RM_WPE_Zon=RM_WPE_Zon;
RM_WPE_Lon=RM_WPE_Lon;
RM_WPE_Mon=RM_WPE_Mon;
RM_WPE_Hon=RM_WPE_Hon;

f6=figure(6);clf
hold all
set(f6, 'Position', [100, 100, 400, 700]);
set(f6, 'DefaultAxesFontSize', 14);
set(f6, 'DefaultLineLineWidth',1.5);
set(f6, 'color','white')    

w=5;
t=1;    

%Process noise
subplot(2,1,1)
hz4=scatter(RM_WPE_Zpn(RM_WPE_Zpn(:,2)==t,w),RMSE_Zpn,30,RM_WPE_Zpn(RM_WPE_Zpn(:,2)==t,1),'o','filled');
        colormap; hold on;
hl4=scatter(RM_WPE_Lpn(RM_WPE_Lpn(:,2)==t,w),RMSE_Lpn,'MarkerEdgeColor',[.6 .6 .6]); hold on; 
hm4=scatter(RM_WPE_Mpn(RM_WPE_Mpn(:,2)==t,w),RMSE_Mpn,'MarkerEdgeColor',[.3 .3 .3]); hold on; 
hh4=scatter(RM_WPE_Hpn(RM_WPE_Hpn(:,2)==t,w),RMSE_Hpn,'MarkerEdgeColor',[0 0 0]); hold on; 
        
ylabel('nRMSE')         
xlabel('WPE')
text(0,5200,'A) Process noise','FontSize', 14, 'FontWeight', 'Bold');        

legend([hh4,hm4,hl4], 'SD = 0.01', 'SD = 0.001', 'SD = 0.0001', 'Location', 'NorthWest');
   legend('Boxoff');
% 
% set(gca, ...  'yticklabel',{'0','0.1','0.2'} , ...
%             'YLim', [-0.01, 0.11]    , ... 
%             'YTick', [0, 0.05, 0.1]      , ...
%             'Linewidth', 1, ...
%             'Box','on'     ) 
        
%Observational noise
subplot(2,1,2)
hz4=scatter(RM_WPE_Zon(RM_WPE_Zon(:,2)==t,w),RMSE_Zpn,30,RM_WPE_Zon(RM_WPE_Zon(:,2)==t,1),'o','filled');
        colormap; hold on;
hl4=scatter(RM_WPE_Lon(RM_WPE_Lon(:,2)==t,w),RMSE_Lon,'MarkerEdgeColor',[.6 .6 .6]); hold on; 
hm4=scatter(RM_WPE_Mon(RM_WPE_Mon(:,2)==t,w),RMSE_Mon,'MarkerEdgeColor',[.3 .3 .3]); hold on; 
hh4=scatter(RM_WPE_Hon(RM_WPE_Hon(:,2)==t,w),RMSE_Hon,'MarkerEdgeColor',[0 0 0]); hold on; 
 
ylabel('nRMSE')         
xlabel('WPE')
text(0,6300,'B) Observational noise','FontSize', 14, 'FontWeight', 'Bold');        

% set(gca, ...  'yticklabel',{'0','0.1','0.2'} , ...
%             'YLim', [-0.01, 0.11]    , ... 
%             'YTick', [0, 0.05, 0.1]      , ...
%             'Linewidth', 1, ...
%             'Box','on'     ) 
        
saveas(f6,'Fig_Ricker_WPE_RMSE_noise','jpg')           
       
%% plot FE as a function of growth rate for different levels of process error

YnRMSE = horzcat(RMSE_Zpn, RMSE_Lpn, RMSE_Mpn); %Forcast error on the y axis
ZGR = horzcat(rs',rs',rs'); %Growth rate on the X axis

f7=figure(7);clf
hold all
set(f7, 'Position', [100, 100, 900, 1000]);
set(f7, 'DefaultAxesFontSize', 14);
set(f7, 'DefaultLineLineWidth',1);
set(f7, 'color','white')   

subplot(3,1,1)
        h3=scatter(bfrctn(1:end-1,1),bfrctn(1:end-1,2)/scale,2,bfrctn(1:end-1,1),'o','filled');
        hold on
       
        set(gca, ... 
            'XLim', [rmin, rmax]         , ...
            'XTick', rmin:0.2:rmax      , ...
            'xticklabel',{'2.4','2.6','2.8','3.0','3.2','3.4'} , ...
            'Linewidth', 1, ...
            'Box','on')
        
         for ii=1:length(r_plot)
            plot([r_plot(ii),r_plot(ii)],[4000,0], '.-', 'Color',arcol(ii,:), 'Linewidth', 2);
        hold on   
         end 
         
        xlabel('Growth rate, r')
        ylabel('Abundance'); 
        text(3.325,1.2,'B)','FontSize', 12, 'FontWeight', 'Bold');
        
subplot(3,1,2)
plot(ZGR, YnRMSE, '-');
ylabel('nRMSE')         
xlabel('Growth rate, r')

subplot(3,1,3)
plot(ZGR, YnRMSE, '-');hold on
plot(rs, RMSE_Hpn,'k-');
ylabel('nRMSE')         
xlabel('Growth rate, r')

legend('No process error', 'Low; SD=0.0001','Med; SD=0.001','High; SD=0.01','Location','northwest')