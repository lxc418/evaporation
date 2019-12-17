%calculate surface resistance
%medium sand

%variable
%evapo_rate_relative_capillary                % relative conductance of vapor through the soil air interface contributed to by the vaporization of the
%                                             % capillary water in the topmost soil layer
%evapo_rate_relative_vapor                    % relative conductance of vapor through the soil air interface contributed to by vapor generated from the
%                                             % vaporization plane beneath the dry soil layer
%Thickness_surface_edl % thickness of the external diffusive layer by surface resistance
%R_0                   % radius of water saturated pores(bottom of funnel)
%R_m                   % radius that corresponds to psi_m according to the Young?Laplace equation (m)
%R_2                   % radius of water saturated pores(top of funnel)
%R_3                   % radius of cylinder building block
%R_c                   % characteristic radius  
%Psi_m                 % matric potential

clear
%constant
xi               = -1.469e-5; % Young_Laplace equation constant(assuming contact angle is zero)
psi_0_m          = -5e4;      % matric potential that corresponds to zero liquid water saturation
tortuosity_0     = 0.66;      % tortuosity of liquid water in soils when the liquid water saturation is zero
diffusivity_m2Ps = 2.62e-5;   % diffusivity of vapor at 22 centigrade
free_path_gas_m  = 0.6e-7;    % mean free path of gas molecules at 22 centigrade

%experimental conditions
thickness_NSL_m      = 0.05;%thickness of the near surface soil layer(NSL)
thickness_aero_edl_m = 5e-3;%thickness of the external diffusive layer(EDL) by aerodynamics

%parameters about soil
psi_p_m              = -10;%matric potential in the NSL corresponding to the initial liquid water saturation at early stage IV(m)
porosity            = 0.39;
saturation_residual = 0.09;%residual liquid water saturation
n                   = 0.0;%correction function between TSL and NSL
beta               = pi/4;%the characteristic angle of soil particle shape

%fitting parameter for the van Genuchten soil water rete3tion curve
% below are working parameters
av_Pm = 10.5;
nv    = 5.5;
radius_particle_m   = 3.5e-4;%average particle size

%%
psim_m_ay= -[0.0001:0.0001:0.01,0.01:0.001:0.1,0.1:0.01:1,1:1:10,10:10:100,100:100:1000,1000:1000:50000];
% psim_m_ay= -[0.001:0.001:100];

r_m_ay   = xi./psim_m_ay;

[saturation_NSL_ay ,saturation_effective_NSL_ay ] = SWCC_Fayer1995WRR(psim_m_ay, -av_Pm, nv, psi_0_m, saturation_residual);

% figure
% semilogx(-psim_m_ay,saturation_NSL_ay)% test swcc
% hold on

r_0_ay                  = -xi*av_Pm*(-av_Pm*psim_m_ay).^(nv-1).*(  (1+ (-av_Pm*psim_m_ay).^-nv)  .^(1-1/nv) -1);
r_c_m                   = radius_particle_m/tan(beta);
r_2_ay                  = r_0_ay+r_c_m;


thickness_funnel_edl_ay = radius_particle_m./exp(1/r_c_m.*(r_m_ay-r_0_ay));
% plot(r_m_ay,-thickness_funnel_edl_ay);

saturation_effective_TSL_ay    = saturation_effective_NSL_ay.^(1+n);
    
% relative_wetted_surface_ay     = saturation_effective_TSL_ay .* (1-thickness_funnel_edl_ay.*(1-porosity)./radius_particle_m);
% relative_wetted_surface_ay    = r_m_ay.^2./r_2_ay.^2;

% relative_wetted_surface_ay     = saturation_effective_TSL_ay .* (r_m_ay./(r_0_ay+radius_particle_m)).^2;
relative_wetted_surface_ay     = saturation_effective_TSL_ay .* ((r_0_ay+radius_particle_m*log(radius_particle_m./thickness_funnel_edl_ay))./(r_0_ay+radius_particle_m)).^2;
% relative_wetted_surface_ay    = saturation_effective_TSL_ay .* porosity;

% water_content_effective_TSL_ay = relative_wetted_surface_ay.*saturation_effective_TSL_ay;
water_content_NSL_ay           = porosity .* saturation_NSL_ay;

r_3_ay                         = r_m_ay./relative_wetted_surface_ay.^0.5;

evapo_rate_relative_capillary = 1./(1  +  r_3_ay.^2 .* thickness_funnel_edl_ay./r_2_ay.^2/thickness_aero_edl_m + ...
          r_m_ay./2./thickness_aero_edl_m./relative_wetted_surface_ay .* ...  
         (  2*free_path_gas_m./r_m_ay  +  1./ (1+free_path_gas_m./r_m_ay ) - relative_wetted_surface_ay.^0.5)  );

% evapo_rate_relative_capillary = 1./(1  +  thickness_funnel_edl_ay./relative_wetted_surface_ay./thickness_aero_edl_m.*(r_m_ay./(r_0_ay+radius_particle_m)).^2 + ...
%           r_m_ay./2./thickness_aero_edl_m./relative_wetted_surface_ay .* ...  
%          (  2*free_path_gas_m./r_m_ay  +  1./ (1+free_path_gas_m./r_m_ay ) - relative_wetted_surface_ay.^0.5)  );

% k_vr_c_ay=1./(1+r_3_ay.^2.*thickness_funnel_edl_ay./r_2_ay.^2/thickness_aero_edl_m);%test thickness influence
% k_vr_c_ay=1./(1+r_3_ay.^2 .* thickness_funnel_edl_ay./r_2_ay.^2/thickness_aero_edl_m + ...
%     free_path/thickness_aero_edl_m./water_content_effective_TSL_ay);                %when free_path>>r_m

evapo_rate_relative_vapor = (thickness_aero_edl_m+thickness_funnel_edl_ay)*tortuosity_0*porosity*log(-psi_0_m) /...
                             thickness_NSL_m./log(-psim_m_ay-psi_p_m);

evapo_rate_relative       = evapo_rate_relative_capillary.*(1-evapo_rate_relative_vapor)+evapo_rate_relative_vapor;

% surface_resistance_ay_sPm =(thickness_funnel_edl_ay+thickness_aero_edl_m)./(diffusivity_m2Ps*evapo_rate_relative)-thickness_aero_edl_m/diffusivity_m2Ps;
surface_resistance_ay_sPm = thickness_aero_edl_m/diffusivity_m2Ps*(1./evapo_rate_relative - 1);
% surface_resistance_ay_sPm = thickness_aero_edl_m/diffusivity_m2Ps./evapo_rate_relative;

%%
figure

% semilogx(-psim_m_ay,water_content_NSL_ay);
% xlim([0.01 inf])

% plot (r_m_ay, -thickness_funnel_edl_ay);
% subplot(2,2,1)
% semilogy(saturation_NSL_ay,evapo_rate_relative_capillary);
% hold on
% hold on;
% semilogy(saturation_NSL_ay,evapo_rate_relative_vapor);
% hold on
% semilogy(saturation_NSL_ay,evapo_rate_relative);
% 
% subplot(2,2,2)
% plot(saturation_NSL_ay,surface_resistance_ay_sPm);

semilogy(saturation_NSL_ay,surface_resistance_ay_sPm);
hold on
ylim([5e-1 inf])
% % xlim([0 0.4])
% 
% % hold on
% subplot(2,2,3)
% plot(water_content_NSL_ay,thickness_funnel_edl_ay);
% ylim([0 2.1e-4])
% 
% subplot(2,2,4)
% semilogy(water_content_NSL_ay,-psim_m_ay);