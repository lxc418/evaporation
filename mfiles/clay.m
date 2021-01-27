clear
%parameters of medium sand

%constant
xi               = -1.46789e-5; % Young_Laplace equation constant(assuming contact angle is zero)
psi_0_m          = -5e4;      % matric potential that corresponds to zero liquid water saturation
tortuosity_0     = 0.66;      % tortuosity of liquid water in soils when the liquid water saturation is zero
diffusivity_m2Ps = 2.63024e-5;   % diffusivity of vapor at 22 centigrade
free_path_gas_m  = 0.6e-7;    % mean free path of gas molecules at 22 centigrade

%experimental conditions
thickness_NSL_m      = 0.03;%thickness of the near surface soil layer(NSL)
thickness_roughness_m = 0.3e-3;%thickness of the external diffusive layer(EDL) by aerodynamics

%parameters about soil
psi_p_m              = -5000;%matric potential in the NSL corresponding to the initial liquid water saturation at early stage IV(m)
porosity             = 0.55;
saturation_residual  = 0.068;%residual liquid water saturation
n                    = 5;%correction function between TSL and NSL
beta                 = pi/4;%the characteristic angle of soil particle shape

%fitting parameter for the van Genuchten soil water rete3tion curve
% below are working parameters
av_Pm                = 0.015;
nv                   = 1.5;
radius_particle_m    = 1e-6;%average particle size

% save medium_sand.mat
figure
new_r_s_model