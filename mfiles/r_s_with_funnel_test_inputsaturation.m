%calculate surface resistance
%coarse sand
clear
%constant
xi               = -1.469e-5; % Young_Laplace equation constant(assuming contact angle is zero)
psi_0_m          = -100e4;      % matric potential that corresponds to zero liquid water saturation
tortuosity_0     = 0.66;      % tortuosity of liquid water in soils when the liquid water saturation is zero
diffusivity_m2Ps = 2.62e-5;   % diffusivity of vapor at 22 centigrade

%experimental conditions
thickness_NSL_m      = 0.05;%thickness of the near?surface soil layer(NSL)
thickness_aero_edl_m = 1e-3;%thickness of the external diffusive layer(EDL) by aerodynamics

%parameters about soil
psi_p               = -10;%matric potential in the NSL corresponding to the initial liquid water saturation at early stage IV(m)
porosity            = 0.40;
saturation_residual = 0.06;%residual liquid water saturation
n                   = 0.5;%correction function between TSL and NSL

%fitting parameter for the van Genuchten soil water retention curve
% av_Pm = 14.5;
% nv    = 2.68; for sand

% below are working parameters
% av_Pm = 14.5;
% nv    = 2.68;
% radius_particle_m   = 0.7e-3;%average particle size


av_Pm = 3.6;
nv    =1.37;
radius_particle_m   = 0.7e-3/8;%average particle size

%variable
%Thickness_surface_edl % thickness of the external diffusive layer by surface resistance
%K_vr_c                % relative conductance of vapor through the soil?air interface contributed to by the vaporization of the
%                      % capillary water in the topmost soil layer
%K_vr_v                % relative conductance of vapor through the soil?air interface contributed to by vapor generated from the
%                      % vaporization plane beneath the dry soil layer
%R_0                   % radius of water saturated pores(bottom of funnel)
%R_m                   % radius that corresponds to psi_m according to the Young?Laplace equation (m)
%R_2                   % radius of water saturated pores(top of funnel)
%R_3                   % radius of cylinder building block
%Psi_m                 % matric potential

%%
water_content_residual      = porosity*saturation_residual;
saturation_NSL_ay           = 0.06:0.0001:1;
water_content_NSL_ay        = porosity*saturation_NSL_ay;

saturation_effective_NSL_ay =(water_content_NSL_ay-water_content_residual)/(porosity-water_content_residual);


psim_m_ay=-(1/av_Pm)*(saturation_effective_NSL_ay.^(nv/(1-nv))-1).^(1/nv);
r_m_ay   = xi./psim_m_ay;


water_content_effective_TSL_ay=porosity*saturation_effective_NSL_ay.^(1+n);


% psim_m_ay= -[0.0001:0.0001:0.01,0.01:0.001:0.1,0.1:0.01:1,1:1:10,10:10:100,100:100:1000,1000:1000:50000];


%opt=SWCC_Fayer1995WRR(psim,-1/o.aa1,o.vn1,-o.phy0,o.swres1);
% [saturation_NSL_ay ,saturation_effective_NSL_ay ] = SWCC_Fayer1995WRR(psim_m_ay, -av_Pm, nv, psi_0_m, saturation_residual);

%figure  
%plot(psim_m_ay,saturation_effective_NSL_ay)
%psim_m_ay = -(saturation_effective_NSL_ay.^(nv/(1-nv))-1).^(1/nv)/av_Pm;%van genuchten

% water_content_effective_TSL_ay = porosity*saturation_effective_NSL_ay.^(1+n);

r_0_ay = -xi*av_Pm*( (1+  (-av_Pm*psim_m_ay).^nv  ).^((nv-1)/nv)-   (-av_Pm*psim_m_ay) .^(nv-1) );%expectation of r from 0 to r_m
r_2_ay = r_0_ay+radius_particle_m;
r_3_ay = r_m_ay./water_content_effective_TSL_ay.^0.5;


r_m_ay(r_m_ay>r_2_ay)  = r_2_ay(r_m_ay>r_2_ay);%if r_m>r_2, r_m=r_2

thickness_funnel_edl_ay=r_2_ay-r_m_ay;

k_vr_c_ay=1./(1+r_3_ay.^2 .* thickness_funnel_edl_ay./r_2_ay.^2/thickness_aero_edl_m  + ...
    r_m_ay/2/thickness_aero_edl_m.* (1./water_content_effective_TSL_ay - 1./water_content_effective_TSL_ay.^0.5));

k_vr_v_ay=(thickness_aero_edl_m+thickness_funnel_edl_ay)*tortuosity_0*porosity*log(-psi_0_m) /...
    thickness_NSL_m./log(-psim_m_ay-psi_p);

k_vr_ay=k_vr_c_ay.*(1-k_vr_v_ay)+k_vr_v_ay;

surface_resistance_ay_sPm=(thickness_funnel_edl_ay+thickness_aero_edl_m)./(diffusivity_m2Ps*k_vr_ay)-thickness_aero_edl_m/diffusivity_m2Ps;
%%
figure
% plot(saturation_NSL_ay,r_m_ay);
% % hold on
% plot(saturation_NSL_ay,water_content_effective_TSL_ay);


subplot(4,1,1)
semilogy(saturation_NSL_ay,k_vr_c_ay);
% hold on;
% semilogy(saturation_NSL_ay,k_vr_v_ay);
% hold on
% semilogy(saturation_NSL_ay,k_vr_ay);

subplot(4,1,2)
semilogy(saturation_NSL_ay,surface_resistance_ay_sPm);

subplot(4,1,3)
plot(saturation_NSL_ay,thickness_funnel_edl_ay);

subplot(4,1,4)
semilogy(saturation_NSL_ay,-psim_m_ay);



