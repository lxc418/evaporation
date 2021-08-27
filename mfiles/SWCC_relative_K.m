%compare swcc from model of Van Genuchten and model of Brooks
clear    
%% Coarse

% AA=14;
% VN=8.5;

% DLAM=5; % Lambda in Corey
% PSI_B=0.06; %sir entry potential

% SWRES=0.06; %residual saturation
% porosity=0.40;
% INTK=10e-9; %intrinsic permeability

%% Medium

% AA=4.5;
% VN=11;
% DLAM=8;
% PSI_B=0.2;
% SWRES=0.09;
% porosity=0.39;
% INTK=5.67e-11;

%% Fine

% AA=4.5;
% VN=8.5;
% DLAM=5.5;
% PSI_B=0.27;
% SWRES=0.1;
% porosity=0.36;
% INTK=1.15e-11;

%% Clay

% AA=0.015;
% VN=1.55;
% 
% DLAM=0.4;
% PSI_B=30;
% 
% PSIR = 436; %https://doi.org/10.1139/cgj-2011-0341
% fred_a = 67.4;
% fred_n = 2.1;
% fred_m = 0.68;
% 
% SWRES=0.1;
% porosity=0.53;
% INTK=1.8e-9; 

%% TEST

AA=14.5;
VN=4.5;
DLAM=1.1;
PSI_B=0.3;
SWRES=0.15;
porosity=0.4;
INTK=1.15e-11;

%%
PSIC0=1000000;
PSIC= [0.01:0.001:1,1:0.01:10,10:1:100,100:10:1000,1000:100:1000000];
      SI       = SWRES*log(PSIC0./PSIC)./log(PSIC0);
      SWRMS1   = 1-SI;
% model of van      
      AAPVN    = 1+(AA*PSIC).^VN;
      VNF      = (VN-1)/VN;                                            
      AAPVNN   = AAPVN.^VNF;
      SE_van   = 1./AAPVNN;
      SW_van   = SI+SWRMS1.*SE_van;
      RELK_van = SE_van.^0.5.*(1-(1-SE_van.^(1/VNF)).^VNF).^2; 
      KK_van = RELK_van*INTK;
      water_content_van = SW_van*porosity;
% model of brooks

      SE_brooks     = (PSI_B./PSIC).^DLAM; %PSI_B air entry value
      SW_brooks     = SWRMS1.*SE_brooks+SI;
      RELK_brooks   = SE_brooks.^(2+2./DLAM+0.5);
      KK_brooks = RELK_brooks*INTK;
      water_content_brooks = SW_brooks*porosity; 
      water_content_brooks(PSIC <= PSI_B) = porosity;
% model of fredlund

%       fred_fitting = 1 - ( log(1+(PSIC./PSIR)) / log(1+(10e5/PSIR)) ); %PSIR suction (kPa) corresponding to residual
%       SW_fred = fred_fitting.*(log (exp(1)+ (PSIC./fred_a).^fred_n) ).^-fred_m;
%       water_content_fred = SW_fred*porosity;


%% plot k.vs.sw      
% figure
% lw=2; %line width
% fz=8; % fontsize
% fl=8; % label font size
% 
% semilogy(SW_brooks, KK_brooks,'LineStyle','--','color',[0.85,0.33,0.10],'LineWidth',lw);
% hold on
% semilogy(SW_van, KK_van,'LineStyle','-','color',[0.85,0.33,0.10],'LineWidth',lw); %et against sat
% 
% 
% xlabel('Saturation(-)','FontSize',fz,'FontWeight','bold')
% ylabel('Hydraulic conductivity of capillary water(-)','FontSize',fz,'FontWeight','bold')
% hleg1 = legend('Coarse sand (Brooks&Corey)','Coarse sand (van Genuchten)','Location','east');
% set(hleg1, 'Box', 'on','FontSize',fz,'FontWeight','bold')
% ax1 = gca;
% set(ax1,'FontSize',fl,'FontWeight','bold')
% pbaspect([1 1 1])
% axis([0 1 10e-35 1])    

%% plot swcc         
figure
lw=2; %line width
fz=8; % fontsize
fl=8; % label font size

semilogx(PSIC,water_content_brooks,'LineStyle','-','color',[0.0 0.45 0.74],'LineWidth',lw); hold on%et against sat
% semilogx(PSIC,SW_brooks,'LineStyle','-.','color',[0.0 0.45 0.74],'LineWidth',lw);hold on %et against sat
% semilogy(SW_brooks,PSIC,'LineStyle','--','color',[0.0 0.45 0.74],'LineWidth',lw); %et against sat
% hold on
semilogy(SW_van,RELK_van,'LineStyle','-.','color',[0.0 0.45 0.74],'LineWidth',lw);

semilogy(PSIC,water_content_van,'LineStyle','-.','color',[0.0 0.45 0.74],'LineWidth',lw);
% semilogx(PSIC,water_content_fred,'LineStyle','--','color',[0.0 0.45 0.74],'LineWidth',lw); %et against sat
% semilogy(SW_van,PSIC,'LineStyle','-','color',[0.0 0.45 0.74],'LineWidth',lw);


xlabel('Matric potential(m)','FontSize',fz,'FontWeight','bold')
ylabel('Volumetric water content(-)','FontSize',fz,'FontWeight','bold')
% hleg1 = legend('Fine sand (Brooks&Corey)','Fine sand (van Genuchten)','Location','Northwest');
% hleg1 = legend('Fine sand (Fredlund)','Fine sand (van Genuchten)','Location','southwest');
% set(hleg1, 'Box', 'on','FontSize',fz,'FontWeight','bold')
ax1 = gca;
set(ax1,'FontSize',fl,'FontWeight','bold')
pbaspect([1 1 1])
% axis([0 10 0 1.05])     


