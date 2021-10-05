ii=1;
% for thickness_diffusion_m = [1e-3 0.1e-3 4e-3 6e-3 10e-3 20e-3]

for thickness_diffusion_m = [ 0.1e-3:0.1e-3:1e-3,1e-3:1e-3:5e-2 ]
    fine_sand
    Rs(ii,:)=surface_resistance_ay;
    MAPE_Rs(ii)=mean(abs((Rs(ii,:)-Rs(10,:))./Rs(ii,:)));
    MaxE_Rs(ii)=max(abs((Rs(ii,:)-Rs(10,:))./Rs(ii,:)));
    MAPE_R(ii)=mean(abs((Rs(ii,:)- Rs(10,:))./(Rs(ii,:)+thickness_diffusion_m/diffusivity_m2Ps )));
    MaxE_R(ii)=max(abs((Rs(ii,:)-Rs(10,:))./(Rs(ii,:)+thickness_diffusion_m/diffusivity_m2Ps )));
    
    ii=ii+1;
end
