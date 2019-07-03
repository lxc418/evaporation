clear
% coarse C9
% calculate k_vr_c
xi=-1.469e-5;
psi_0=-5e4;
delta_a=1e-3;
 
% n=0.5;
theta_p=0.40;
r_p=0.7e-3;
saturation_lr=0.06;
psi_p=-10;
% psi_b=-0.06;
av=14.5;
nv=2.68;

theta_lr=saturation_lr*theta_p;


% M=[];
% r_M=[];
% r_S=[];
% r_S1=[];
% r_S2=[];
% K_vr_cn=[];
Delta_s=[];
K_vr_cf=[];
r_S=[];
R_0=[];
for saturation_l=0.061:0.001:1

    theta=saturation_l*theta_p;

    psi_m=(1/av)*(((theta-theta_lr)/(theta_p-theta_lr))^(nv/(1-nv))-1)^(1/nv);

    saturation_le=(1+(av*psi_m)^nv)^((1-nv)/nv);

    theta_w=theta_p*saturation_le;
    r_0=-xi*av*((1+(av*psi_m)^nv)^((nv-1)/nv)-(av*psi_m)^(nv-1));

R_0=[R_0,r_0];
k_vr_cf=1/(1+(r_0/delta_a+r_p/delta_a)*(1/2/theta_w^0.5-theta_w+0.5));%r_p=delta_0

K_vr_cf=[K_vr_cf,k_vr_cf];
delta_s=theta_w^0.5*(r_0+r_p)*(1/theta_w^0.5-1);
Delta_s=[Delta_s,delta_s];

%
tau_0=0.66;
l_0=0.05;


k_vr_v=((delta_a+delta_s)*tau_0*theta_p*log(-psi_0))/(l_0*log(-psi_m-psi_p));

k_vr=k_vr_cf*(1-k_vr_v)+k_vr_v;

% surface resistance
d_v=2.62e-5; %22.c
r_s=(delta_s)/(d_v*k_vr);

r_S=[r_S,r_s];

end

% n=0.1:0.05:1;
% % semilogy (n,r_S);
% % hold on;
% % semilogy (n,r_S1);
% % hold on;
% % semilogy (n,r_S2);
% % 
% % n=0.1:0.01:1;
% % semilogy (n,KK);
% 
% 
% % semilogy (n,K_vr_cn);
% hold on;
% semilogy (n,K_vr_cf);
% 
% 
% xlim([0.09 1.01]);
% ylim([0.5 1.01]);
% set(gca,'XDir','reverse');   