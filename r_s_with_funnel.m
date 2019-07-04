%calculate surface resistance
clear
%constant
xi=-1.469e-5;% Young_Laplace equation constant(assuming contact angle is zero)
psi_0=-5e4;%matric potential that corresponds to zero liquid water saturation
tortuosity_0=0.66;%tortuosity of liquid water in soils when the liquid water saturation is zero
diffusivity=2.62e-5; %diffusivity of vapor at 22 centigrade

%experimental conditions
thickness_nsl=0.05;%thickness of the near?surface soil layer(NSL)
thickness_aero_edl=1e-3;%thickness of the external diffusive layer(EDL) by aerodynamics

%parameters about soil
radius_particle=0.7e-3;%average particle size
psi_p=-10;%matric potential in the NSL corresponding to the initial liquid water saturation at early stage IV(m)
porosity=0.40;
saturation_residual=0.06;%residual liquid water saturation

%fitting parameter for the van Genuchten soil water retention curve
av=14.5;
nv=2.68;

n=0.5;%correction function between TSL and NSL

%variable
Thickness_surface_edl=[];%thickness of the external diffusive layer by surface resistance
K_vr_c=[];%relative conductance of vapor through the soil?air interface contributed to by the vaporization of the 
           %capillary water in the topmost soil layer
K_vr_v=[];%relative conductance of vapor through the soil?air interface contributed to by vapor generated from the 
           %vaporization plane beneath the dry soil layer
Surface_resistance=[];
R_0=[];%radius of water saturated pores(bottom of funnel)
R_m=[];%radius that corresponds to psi_m according to the Young?Laplace equation (m)
R_2=[];%radius of water saturated pores(top of funnel)
R_3=[];%radius of cylinder building block
Psi_m=[];%matric potential

for saturation=0.06:0.001:1%from residual saturation to 1

    water_content_residual=saturation_residual*porosity;
    water_content=saturation*porosity;
    saturation_effective=(water_content-water_content_residual)/(porosity-water_content_residual);
    
    psi_m=-1/av*((saturation_effective)^(nv/(1-nv))-1)^(1/nv);%van genuchten
    r_m=xi/psi_m;
    
Psi_m=[Psi_m,psi_m];
R_m=[R_m,r_m];    

    water_content_effective=porosity*saturation_effective^(1+n);
    
    r_0=-xi*av*((1+(av*-psi_m)^nv)^((nv-1)/nv)-(av*-psi_m)^(nv-1));%expectation of r from 0 to r_m
    r_2=r_0+radius_particle;
    r_3=r_m/water_content_effective^0.5;
    
R_0=[R_0,r_0];
R_2=[R_2,r_2];
R_3=[R_3,r_3];

if r_m>r_2
   r_m=r_2;
end
   thickness_surface_edl=r_2-r_m;
Thickness_surface_edl=[Thickness_surface_edl,thickness_surface_edl];

   k_vr_c=1/(1+(r_3^2*thickness_surface_edl/r_2^2/thickness_aero_edl)+...
       (r_m/2/thickness_aero_edl)*(1/water_content_effective-1/water_content_effective^0.5));
   k_vr_v=((thickness_aero_edl+thickness_surface_edl)*tortuosity_0*porosity*log(-psi_0))/...
       (thickness_nsl*log(-psi_m-psi_p));
   k_vr=k_vr_c*(1-k_vr_v)+k_vr_v;

K_vr_c=[K_vr_c,k_vr_c];
K_vr_v=[K_vr_v,k_vr_v];

   surface_resistance=(thickness_surface_edl+thickness_aero_edl)/(diffusivity*k_vr)-thickness_aero_edl/diffusivity;

Surface_resistance=[Surface_resistance,surface_resistance];

end

%plot
ii=0.06:0.001:1;
semilogy(ii,Surface_resistance);
hold on;
