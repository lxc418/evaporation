%calculate surface resistance
%coarse sand
clear
%constant
xi           = -1.469e-5; % Young_Laplace equation constant(assuming contact angle is zero)
psi_0_m        = -5e4;      % matric potential that corresponds to zero liquid water saturation
tortuosity_0 = 0.66;      % tortuosity of liquid water in soils when the liquid water saturation is zero
diffusivity_m2Ps  = 2.62e-5;   % diffusivity of vapor at 22 centigrade

%experimental conditions
thickness_nsl      = 0.05;%thickness of the near?surface soil layer(NSL)
thickness_aero_edl = 1e-3;%thickness of the external diffusive layer(EDL) by aerodynamics

%parameters about soil
radius_particle_m   = 0.7e-3;%average particle size
psi_p               = -10;%matric potential in the NSL corresponding to the initial liquid water saturation at early stage IV(m)
porosity            = 0.40;
saturation_residual = 0.06;%residual liquid water saturation

%fitting parameter for the van Genuchten soil water retention curve
av = 14.5;
nv = 2.68;

n=0.5;%correction function between TSL and NSL

%variable
Thickness_surface_edl=[]; % thickness of the external diffusive layer by surface resistance
K_vr_c=[];                % relative conductance of vapor through the soil?air interface contributed to by the vaporization of the
                          % capillary water in the topmost soil layer
K_vr_v=[];                % relative conductance of vapor through the soil?air interface contributed to by vapor generated from the
                          % vaporization plane beneath the dry soil layer
%Surface_resistance = [];
%R_0                = []; % radius of water saturated pores(bottom of funnel)
%R_m                = []; % radius that corresponds to psi_m according to the Young?Laplace equation (m)
%R_2                = []; % radius of water saturated pores(top of funnel)
%R_3                = []; % radius of cylinder building block
%Psi_m              = []; % matric potential

saturation_ay           = 0.06:0.001:1;
water_content_residual  = saturation_residual*porosity;
water_content_ay        = saturation_ay*porosity;
saturation_effective_ay = (water_content_ay-water_content_residual)/(porosity-water_content_residual);


%figure  
%plot(psi_m_ay,saturation_effective_ay)


psi_m_ay = -1/av*(saturation_effective_ay.^(nv/(1-nv))-1).^(1/nv);%van genuchten
r_m_ay   = xi./psi_m_ay;

water_content_effective_NSL_ay = porosity*saturation_effective_ay.^(1+n);

r_0_ay = -xi*av*( (1+  (-av*psi_m_ay).^nv  ).^((nv-1)/nv)-   (-av*psi_m_ay) .^(nv-1) );%expectation of r from 0 to r_m
r_2_ay = r_0_ay+radius_particle_m;
r_3_ay = r_m_ay./water_content_effective_NSL_ay.^0.5;


r_m_ay(r_m_ay>r_2_ay)  = r_2_ay(r_m_ay>r_2_ay);

thickness_surface_edl_ay=r_2_ay-r_m_ay;


k_vr_c_ay=1./(1+r_3_ay.^2 .* thickness_surface_edl_ay./r_2_ay.^2/thickness_aero_edl  + ...
    r_m_ay/2/thickness_aero_edl.* (1./water_content_effective_NSL_ay-1./water_content_effective_NSL_ay.^0.5));


k_vr_v_ay=(thickness_aero_edl+thickness_surface_edl_ay)*tortuosity_0*porosity*log(-psi_0_m) /...
    thickness_nsl./log(-psi_m_ay-psi_p);

k_vr_ay=k_vr_c_ay.*(1-k_vr_v_ay)+k_vr_v_ay;


surface_resistance_ay_sPm=(thickness_surface_edl_ay+thickness_aero_edl)./(diffusivity_m2Ps*k_vr_ay)-thickness_aero_edl/diffusivity_m2Ps;

semilogy(saturation_ay,surface_resistance_ay_sPm);
hold on;

