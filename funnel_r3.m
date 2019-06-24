clear
% fine sandy F10
% calculate k_vr_c

delta_a=2e-3;
r_0=1e-4;


K_vr_1=[];

for theta_1=0:0.001:1
    
    r_p=0.7e-3;
    
    k_vr_1=1/(1+(r_0/delta_a+r_p/delta_a)*(1/(2*sqrt(theta_1))-sqrt(theta_1)+0.5));
    K_vr_1=[K_vr_1,k_vr_1];

end

K_vr_2=[];
for theta_2=0:0.001:1
    
    r_p=0.35e-3;
    
    k_vr_2=1/(1+(r_0/delta_a+r_p/delta_a)*(1/(2*sqrt(theta_2))-sqrt(theta_2)+0.5));
    K_vr_2=[K_vr_2,k_vr_2];

end

K_vr_3=[];
for theta_3=0:0.001:1
    
    r_p=0.2e-4;
    
    k_vr_3=1/(1+(r_0/delta_a+r_p/delta_a)*(1/(2*sqrt(theta_3))-sqrt(theta_3)+0.5));
    K_vr_3=[K_vr_3,k_vr_3];

end


n=0:0.001:1;

semilogy (n,K_vr_1);
hold on;
semilogy (n,K_vr_2);
hold on;
semilogy (n,K_vr_3);
legend('coarse','medium','fine');

xlim([0.09 1.01]);
ylim([0.5 1.01]);
set(gca,'XDir','reverse');      