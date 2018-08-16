function [T_vec, K_result, Ca_result, Na_result, Cl_result, glu_result, gaba_result] = genTuckwell2D(ratio)
dt = 5e-3;
dx = 1e-2;
max_t = 6;
L = 1;
[xv, yv] = meshgrid(-L:dx:L, -L:dx:L);
K_oR = 3.;
Ca_oR = 1.;
Na_oR = 120.;
Cl_oR = 136.25;

u1 = (ones(size(xv)) * K_oR);
u2 = (ones(size(xv)) * Ca_oR);
u3 = (ones(size(xv)) * Na_oR);
u4 = (ones(size(xv)) * Cl_oR);
u5 = (zeros(size(xv)));
u6 = u5;
patch = 17 * exp(-(xv.^2 + yv.^2)/ (0.05^2));
u1=u1+patch;
u4=u4+patch;

T_vec=0:dt:max_t;
u0 = cat(3, u1, u2, u3, u4, u5, u6);
u=u0;
K_result=(zeros([size(xv),length(T_vec)]));
K_result(:,:,1)=(u1);
Ca_result=(zeros(size(K_result)));
Ca_result(:,:,1)=(u2);
Na_result=(zeros(size(K_result)));
Na_result(:,:,1)=(u3);
Cl_result=(zeros(size(K_result)));
Cl_result(:,:,1)=(u4);
glu_result=(zeros(size(K_result)));
glu_result(:,:,1)=(u5);
gaba_result=(zeros(size(K_result)));
gaba_result(:,:,1)=(u6);

% Define constants
% Ratios of extracellular to intracellular
alpha1 = 0.25;
alpha2 = 2.0;

% Diffusion coefficients (cm2 s-1)
D_K = 2.5e-3;
D_Ca = 1.e-3;
D_Na = 1.7e-3;
D_Cl = 2.5e-3;
D_Te = 1.3e-3;
D_Ti = 1.3e-3;

% Resting concentration (mM)
K_oR = 3.;
K_iR = 140.;
Ca_oR = 1.;
Ca_iR = 1e-4;
Na_oR = 120.;
Na_iR = 15.;
Cl_oR = 136.25;
Cl_iR = 6.;

% Permeabilities
p_Na = 0.05;
p_Cl = 0.4;

% Ca conductance params (mV)
VM_star = 45.;
VM_T = -60.;

% Dynamics
% k1 = 78.091;
k1 = 90;
k2 = 1.5;
k3 = 0.;
k4 = 1.5;
k5 = 108.2604;
k6 = 1.5e-4;
k7 = 0.2;
k8 = 3.998e-4;
k9 = 2;
k10 = 0.;
k11 = 39.814;
k12 = -104.05;
k13 = 0.;
k14 = 104.064;
k15 = -3.47;
k16 = -3.15;
k17 = 577.895*ratio;
k18 = 2.5;
k19 = 2.5;
k20 = 0.8;
k21 = 0.2;
k22 = 0.3677;
k25 = 260.16;
k26 = 9.0;
k27 = 47.124*ratio;
k28 = 1.;
k29 = 47.124*ratio;
k30 = 1.;
k31 = 0.11;

% Misc.
F = 9.649e4;  % C mol-1
R = 8.31;  % V C mol-1 K-1
T = 310.;  % K

cnt = 2;
for time = dt:dt:max_t
    K_o = u(:,:,1);
    Ca_o = u(:,:,2);
    Na_o = u(:,:,3);
    Cl_o = u(:,:,4);
    T_E = u(:,:,5);
    T_I = u(:,:,6);

    % The model
    K_i = K_iR + alpha1 * (u1 - K_o);
    Ca_i = Ca_iR + alpha2 * (u2 - Ca_o);
    Na_i = Na_iR + alpha1 * (u3 - Na_o);
    Cl_i = Cl_iR + alpha1 * (u4 - Cl_o);

    Vm = 1000*log((K_o + p_Na * Na_o + p_Cl * Cl_i)./(K_i + p_Na * Na_i + p_Cl * Cl_o)) * R * T / F;
    VmR = 1000*log((K_oR + p_Na * Na_oR + p_Cl * Cl_iR)./(K_iR + p_Na * Na_iR + p_Cl * Cl_oR)) * R * T / F;
    Vk = 1000 * log(K_o./K_i) * R * T / F;
    Vca = 1000 * 0.5 * log(Ca_o./Ca_i) * R * T / F;
    Vna = 1000 * log(Na_o./Na_i) * R * T / F;
    Vcl = -1000 * log(Cl_o ./ Cl_i) * R * T / F;

    f_Kp = k6 * (Vm - VmR) .* (Vm - Vk) .* heaviside(Vm - VmR);
    P_KNa = k17 * (1 + k18./Na_i).^(-3) .* (1 + k19./K_o).^(-2).*heaviside(Na_i).*heaviside(K_o);

    k5 = k17 * (1 + k18./Na_iR).^(-3) .* (1 + k19./K_oR).^(-2);
    f_K = k1 * (Vm - Vk) .* (T_E.*heaviside(T_E)./(T_E+k2) + k3*T_I.*heaviside(T_I)./(T_I+k4)) - P_KNa +f_Kp + k5;    

    k32 = 1 + tanh(k31 * (VM_T + VM_star));
    g_ca = (1 + tanh(k31 * (Vm + VM_star)) - k32) .* heaviside(Vm - VM_T);
    P_Ca = k20 * Ca_i .* heaviside(Ca_i) ./ (Ca_i + k21);
    k8 = k20 * Ca_iR ./ (Ca_iR + k21);
    f_Ca = k7 * (Vm - Vca) .* g_ca + P_Ca -k8;

    k11 = k22 * k5;
    f_Na = k9 * (Vm - Vna) .* (T_E.*heaviside(T_E)./(T_E+k2) + k10*T_I.*heaviside(T_I)./(T_I+k4)) + k22 * P_KNa - k11;

    P_Cl = k25 * Cl_i .* heaviside(Cl_i) ./ (Cl_i + k26);
    k14 = k25 * Cl_iR ./ (Cl_iR + k26);
    f_Cl = k12 * (Vm - Vcl) .* (k13*T_E.*heaviside(T_E)./(T_E+k2) + T_I.*heaviside(T_I)./(T_I+k4)) + P_Cl - k14;

    f_TE = k15 * (Vm - Vca) .* g_ca - k27 * T_E .* heaviside(T_E) ./ (T_E + k28);
    f_TI = k16 * (Vm - Vca) .* g_ca - k29 * T_I .* heaviside(T_I) ./ (T_I + k30);

    beta1=dt*D_K/(dx*dx);
    beta2=dt*D_Ca/(dx*dx);
    beta3=dt*D_Na/(dx*dx);
    beta4=dt*D_Cl/(dx*dx);
    beta5=dt*D_Te/(dx*dx);
    beta6=dt*D_Ti/(dx*dx);

    u(2:end-1,2:end-1,1)=u(2:end-1,2:end-1,1)+beta1*(u(1:end-2,2:end-1,1)+u(3:end,2:end-1,1)+...
        u(2:end-1,1:end-2,1)+u(2:end-1,3:end,1)-4*u(2:end-1,2:end-1,1))+dt*f_K(2:end-1,2:end-1,1);
    u(1,:,1)=K_oR;
    u(end,:,1)=K_oR;
    u(:,1,1)=K_oR;
    u(:,end,1)=K_oR;

    u(2:end-1,2:end-1,2)=u(2:end-1,2:end-1,2)+beta2*(u(1:end-2,2:end-1,2)+u(3:end,2:end-1,2)+...
        u(2:end-1,1:end-2,2)+u(2:end-1,3:end,2)-4*u(2:end-1,2:end-1,2))+dt*f_Ca(2:end-1,2:end-1,1);
    u(1,:,2)=Ca_oR;
    u(end,:,2)=Ca_oR;
    u(:,1,2)=Ca_oR;
    u(:,end,2)=Ca_oR;

    u(2:end-1,2:end-1,3)=u(2:end-1,2:end-1,3)+beta3*(u(1:end-2,2:end-1,3)+u(3:end,2:end-1,3)+...
        u(2:end-1,1:end-2,3)+u(2:end-1,3:end,3)-4*u(2:end-1,2:end-1,3))+dt*f_Na(2:end-1,2:end-1,1);
    u(1,:,3)=Na_oR;
    u(end,:,3)=Na_oR;
    u(:,1,3)=Na_oR;
    u(:,end,3)=Na_oR;

    u(2:end-1,2:end-1,4)=u(2:end-1,2:end-1,4)+beta4*(u(1:end-2,2:end-1,4)+u(3:end,2:end-1,4)+...
        u(2:end-1,1:end-2,4)+u(2:end-1,3:end,4)-4*u(2:end-1,2:end-1,4))+dt*f_Cl(2:end-1,2:end-1,1);
    u(1,:,4)=Cl_oR;
    u(end,:,4)=Cl_oR;
    u(:,1,4)=Cl_oR;
    u(:,end,4)=Cl_oR;

    u(2:end-1,2:end-1,5)=u(2:end-1,2:end-1,5)+beta5*(u(1:end-2,2:end-1,5)+u(3:end,2:end-1,5)+...
        u(2:end-1,1:end-2,5)+u(2:end-1,3:end,5)-4*u(2:end-1,2:end-1,5))+dt*f_TE(2:end-1,2:end-1,1);
    u(1,:,5)=0;
    u(end,:,5)=0;
    u(:,1,5)=0;
    u(:,end,5)=0;

    u(2:end-1,2:end-1,6)=u(2:end-1,2:end-1,6)+beta6*(u(1:end-2,2:end-1,6)+u(3:end,2:end-1,6)+...
        u(2:end-1,1:end-2,6)+u(2:end-1,3:end,6)-4*u(2:end-1,2:end-1,6))+dt*f_TI(2:end-1,2:end-1,1);
    u(1,:,6)=0;
    u(end,:,6)=0;
    u(:,1,6)=0;
    u(:,end,6)=0;

    if ~isreal(u)
        error('result becomes imaginary!')
    end
    
    K_result(:,:,cnt)=(u(:,:,1));
    Ca_result(:,:,cnt)=(u(:,:,2));
    Na_result(:,:,cnt)=(u(:,:,3));
    Cl_result(:,:,cnt)=(u(:,:,4));
    glu_result(:,:,cnt)=(u(:,:,5));
    gaba_result(:,:,cnt)=(u(:,:,6));
    
    cnt=cnt+1;
end

figure,plot(squeeze(K_result(101,151,:)));
end

function out=heaviside(in)
out = double(in>0);
end
