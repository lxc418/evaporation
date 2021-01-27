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

%%
%psim_m_ay= -[0.0001:0.0001:0.01,0.01:0.001:0.1,0.1:0.001:1,1:0.1:10,10:1:100,100:10:1000,1000:500:50000];
psim_m_ay= -[0.001:0.001:1,1:0.1:10,10:1:100,100:10:1000,1000:500:100000];

r_m_ay   = xi./psim_m_ay;
r_0_ay   = -xi*av_Pm*(-av_Pm*psim_m_ay).^(nv-1).*(  (1+ (-av_Pm*psim_m_ay).^-nv)  .^(1-1/nv) -1);
r_c_m    = radius_particle_m/tan(beta);
r_2_ay   = r_0_ay+r_c_m;
thickness_funnel_ay  = radius_particle_m./exp(1/r_c_m.*(r_m_ay-r_0_ay));

%reference evaporation
% psim_m_ref = -0.1;
% r_m_ref    = xi/psim_m_ref;
% r_0_ref    = -xi*av_Pm*(-av_Pm*psim_m_ref).^(nv-1).*(  (1+ (-av_Pm*psim_m_ref).^-nv)  .^(1-1/nv) -1);
% r_2_ref    = r_0_ref+r_c_m;
% thickness_funnel_ref  = radius_particle_m./exp(1/r_c_m.*(r_m_ref-r_0_ref));


[saturation_NSL_ay ,saturation_effective_NSL_ay ] = SWCC_Fayer1995WRR(psim_m_ay, -av_Pm, nv, psi_0_m, saturation_residual);
% figure
% semilogx(-psim_m_ay,saturation_NSL_ay)% test swcc
% hold on

% plot(r_m_ay,-thickness_funnel_ay);
saturation_effective_TSL_ay    = saturation_effective_NSL_ay.^(1+n);
r_m_ay(r_m_ay>r_2_ay)          = r_2_ay(r_m_ay>r_2_ay);%if r_m>r_2, r_m=r_2    
relative_wetted_surface_ay     = saturation_effective_TSL_ay .* (1-thickness_funnel_ay.*(1-porosity)./radius_particle_m);
% relative_wetted_surface_ay     = saturation_effective_TSL_ay .* (r_m_ay.^2./r_2_ay.^2);
% saturation_TSL_ay              = saturation_NSL_ay.^(1+n);
% water_content_TSL_ay           = porosity .* saturation_TSL_ay;
water_content_NSL_ay           = porosity .* saturation_NSL_ay;
r_3_ay                         = r_m_ay./(relative_wetted_surface_ay.^0.5);

r_3_ay(r_3_ay<r_2_ay)          = r_2_ay(r_3_ay<r_2_ay);
r_m_ay(r_3_ay<r_m_ay)          = r_3_ay(r_3_ay<r_m_ay);
evapo_rate_relative_capillary  = 1./(1  +  r_3_ay.^2 .* thickness_funnel_ay./r_2_ay.^2/thickness_roughness_m + ...
                                r_m_ay./2./thickness_roughness_m./relative_wetted_surface_ay .* ...  
                               ( 2*free_path_gas_m./r_m_ay  +  1./ (1+free_path_gas_m./r_m_ay ) - relative_wetted_surface_ay.^0.5)  );

% evapo_rate_relative_capillary  = ((1./r_m_ref-1./r_2_ref+ 2*thickness_funnel_ref/r_2_ref.^2) )  .*   (1./( (r_3_ay-r_m_ay)./(r_m_ay.*r_3_ay) + ...
%                                   2*thickness_funnel_ay./r_2_ay.^2));
     
thickness_dry_ay               = thickness_NSL_m*log(-psim_m_ay-psi_p_m)./log(-psi_0_m);
                           
evapo_rate_relative_vapor      = 1./(thickness_dry_ay/(tortuosity_0*porosity*thickness_roughness_m)+1);
% evapo_rate_relative_vapor      = (tortuosity_0*porosity*thickness_roughness_m)./thickness_dry_ay;

% evapo_rate_relative_vapor      = ((1./r_m_ref-1./r_2_ref+ 2*thickness_funnel_ref/r_2_ref.^2) )    .*   ( 1./((r_3_ay-r_0_ay)./(r_0_ay.*r_3_ay) + ...
%                                   2*thickness_dry_ay./r_3_ay.^2./tortuosity_0));
                    

evapo_rate_relative            = evapo_rate_relative_capillary.*(1-evapo_rate_relative_vapor)+evapo_rate_relative_vapor;

% surface_resistance_ay          = 1./evapo_rate_relative;
% surface_resistance_ay          = thickness_roughness_m/diffusivity_m2Ps*(1./evapo_rate_relative - 1);
surface_resistance_ay          = thickness_roughness_m/diffusivity_m2Ps./evapo_rate_relative-thickness_roughness_m/diffusivity_m2Ps;

surface_resistance_van         = 10*exp(35.63*(0.15-water_content_NSL_ay));

evapo_rate_relative_Sch        = 1./(1+2/pi*radius_particle_m/thickness_roughness_m.*(pi/4./water_content_NSL_ay).^0.5.*((pi/4./water_content_NSL_ay).^0.5-1));
surface_resistance_Sch         = thickness_roughness_m/diffusivity_m2Ps./evapo_rate_relative_Sch;
% surface_resistance_sun         = 0.335+0.035*(porosity./water_content_NSL_ay).^2.3;
%%
% figure

% semilogx(-psim_m_ay,water_content_NSL_ay);
% xlim([0.01 inf])

% plot (r_m_ay, -thickness_funnel_ay);
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

semilogy(saturation_NSL_ay,surface_resistance_ay);
hold on
% semilogy(saturation_NSL_ay,surface_resistance_van);
% semilogy(saturation_NSL_ay,surface_resistance_Sch);
% legend ('new', 'van', 'Sch')
% % semilogy(saturation_NSL_ay,surface_resistance_sun);

% ylim([5e-1 10e5])
% % xlim([0 0.4])
% 
% % hold on
% subplot(2,2,3)
% plot(water_content_NSL_ay,thickness_funnel_ay);
% ylim([0 2.1e-4])
% 
% subplot(2,2,4)
% semilogy(water_content_NSL_ay,-psim_m_ay);