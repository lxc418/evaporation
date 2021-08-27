%% Chenming's surface resistance 
psim_m_cm= -[-psi_b:0.001:1,1:0.1:10,10:1:100,100:10:1000,1000:500:-psi_0_m];

thickness_dry_ay               = thickness_NSL_m*log(-psim_m_cm-psi_p_m)./log(-psi_0_m);
evapo_rate_relative_vapor      = 1./(thickness_dry_ay/(tortuosity_0*porosity*thickness_diffusion_m)+1);

saturation_effective_NSL_cm=(psi_b./psim_m_cm).^DLAM;
saturation_effective_TSL_cm=saturation_effective_NSL_cm.^(1+n);
beta = (log(-psi_0_m)-log(-psim_m_cm))/log(-psi_0_m);
saturation_NSL_cm = (1-beta*saturation_residual).*saturation_effective_NSL_cm +beta*saturation_residual;

PORSE=porosity*saturation_effective_TSL_cm;
%      PARAMETER IN THE HYPERGEOMETRIC
C =(1./PORSE-1./sqrt(PORSE))./2./thickness_diffusion_m;
R1EFF=DLAM*xi/(DLAM+1)./psim_m_cm;
RELKV=1./(1+C.*R1EFF).*(1-evapo_rate_relative_vapor)+evapo_rate_relative_vapor;
surface_resistance_cm=thickness_diffusion_m/diffusivity_m2Ps./RELKV-thickness_diffusion_m/diffusivity_m2Ps;

%% van de Griend & Owe
surface_resistance_van         = 10*exp(35.63*(0.15-water_content_NSL_ay));

%% Schlunder
evapo_rate_relative_Sch        = 1./(1+2/pi*radius_particle_m/thickness_diffusion_m.*(pi/4./water_content_NSL_ay).^0.5.*((pi/4./water_content_NSL_ay).^0.5-1));
% surface_resistance_Sch         = thickness_diffusion_m/diffusivity_m2Ps./evapo_rate_relative_Sch;
surface_resistance_Sch          = thickness_diffusion_m/diffusivity_m2Ps./evapo_rate_relative_Sch-thickness_diffusion_m/diffusivity_m2Ps;

%% Shu Fen Sun
% surface_resistance_sun         = 0.335+0.035*(porosity./water_content_NSL_ay).^2.3;
