clear;
%% Tuckwell model
percent = 100;  % Reduce efficiency of K and neurotransmitter cleaning up
[T_vec, K_result, Ca_result, Na_result, Cl_result, glu_result, gaba_result] = genTuckwell2D(percent/100);
save(['tuckwell_',num2str(percent),'percent'], 'Ca_result', 'Cl_result', 'K_result', 'Na_result', 'gaba_result', 'glu_result', 'T_vec');

%% Neurovascular coupling model
% First add some baseline before KCl application
dt = T_vec(2) - T_vec(1);
T = 26*[0:dt:3.5-dt, 3.5+T_vec]'; % Scale time axis
K = cat(3, K_result(1)*ones([201,201,3.5/dt]), K_result);
Glu = cat(3, glu_result(1)*ones([201,201,3.5/dt]), glu_result);
clear T_* K_* glu_*
% Take a slice of the pizza
tmp = blkdiag(tril(ones(101)), zeros(100));
K_slice = zeros(length(T), 5151);
Glu_slice = K_slice;
cnt = 1;
for i=1:201*101
    if tmp(i)>0
        [x,y] = ind2sub([201,201],i);
        K_slice(:,cnt) = K(x,y,:);
        Glu_slice(:,cnt) = Glu(x,y,:);
        cnt=cnt+1;
    end
end
save(['input_patch_', num2str(percent)],'T','K_slice', 'Glu_slice');
% Take a bite
K_bite = squeeze(K(101, 131, :));
Glu_bite = squeeze(Glu(101,131,:));
% Set up the NV coupling model
Ca_astr = 150.;          % nM
Ca_ER = 4.e5;             % nM
IP3 = 0.01 ;             % uM
EET = 0.4 ;                % uM
h = 0.1/(0.1+Ca_astr);               % NO UNIT
Vk = -85.0 ;               % mV
v4 = 14.5;           % mV
v5 = 8.0;            % mV
v6 = -15.0 ;         % mV
Ca3 = 400.;           % uM
Ca4 = 150.;          % uM
phi_n = 2.664 ;      % s-1

Vm=-45.;             % mV
Ca=500.;             % nM

v3=-0.5*v5*tanh((Ca-Ca3)/Ca4)+v6;
n=0.5*(1+tanh((Vm-v3)/v4));

v3 = -0.5 * v5 * tanh((Ca_astr - Ca3) / Ca4) + v6;
lambda_BK = phi_n * cosh((Vk - v3) / 2 / v4);
n_BK = 0.5 * (1 + tanh((Vk + 2. * EET - v3) / v4));

t = (0:0.05:245)'; % The new time axis
INIT = [Ca_astr, Ca_ER, IP3, h, EET, n_BK, 100., Vk, n, Vm, Ca, 0., 5., 0.5, 0.7, 0.1709];

opt=odeset('MaxStep',0.1,'BDF','on');
[T,result]=ode15s(@(t,y)nvcoupling(t, y, K_bite, Glu_bite, T, 80, 1),t,INIT,opt); % CPP=80mmHg
R = result(:,end)/2/pi;
figure,plot(t,R)
Kp = result(:, 13);
save(['ProbedRadius_101_131_', num2str(percent)], 'T', 'R', 'Kp');
% clear;
% load(['ProbedRadius_101_131_', num2str(percent)]);

%% And hemodynamics
T0=100; %%% uM
Tc      = 1.2;  %%% capillary transit time
Tv      = 6;    %%% venous transit time
pa      = 0.1;        %%% CBV_a_0 / CBV_0 ... arterial contribution. 0.04 means 4%
pc      = 0.4;  %%% F*CBV_c_0 / CBV_0 ... capillary contribution with Fahreaus factor 
alpha   = 0.56; %%% rate constant 
pv = 1-pa-pc;

Sa    = 0.98;          %arterial saturation
Far   = 0.8;           %Fahraeus factor

Sc    = Sa*(1/Tc/alpha)*(1-exp(-alpha*Tc));            %average capillary saturation
Sv    = Sa*exp(-alpha*Tc);            %venous saturation
ctHb=2300;

phi=T0/ctHb; 
phi_a  = pa*phi; 
phi_cF = pc*phi; 
phi_c=phi_cF/Far;
phi_v  = (1-(pa+pc))*phi;

R = smooth(R, 100); % smooth with a 5sec window!
R0 = mean(R(801:1600));    % cm

va_t = (R/R0).^2 - 1;
va_t=va_t';

vv_t = va_t;
cbv = vv_t * (phi_a/phi+phi_v/phi); 

fs=1/(T(2)-T(1)); % Hz, sample rate

cbf = (cbv + 1).^(1/0.38)-1;

% Get metabolism
E0 = 0.7;
cmro2 = (1-(1-E0).^(1./(1+cbf)))/E0.*(1+cbf) - 1;

% Sergio's model
b1 = exp(1)/Tc*exp(-exp(1)*T/Tc);
b2 = 1/(0.6*(Tc+Tv))*exp(-pi*(T-0.5*(Tc+Tv)).^2/(0.6*(Tc+Tv))^2);
D0 = ctHb*((1-Sa)*phi_a + (1-Sc)*Far*phi_c + (1-Sv)*phi_v);
O0 = ctHb*(Sa*phi_a + Sc*Far*phi_c + Sv*phi_v);
HbD = ctHb*((1-Sa)*phi_a*(1+va_t) + (1-Sc)*Far*phi_c + (1-Sv)*phi_v*(1+vv_t) ...
    - (Sc/Sv*(Sc - Sv)*Far*phi_c*filter(b1,1,cbf - cmro2)/fs + (Sa-Sv)*phi_v*filter(b2,1,cbf - cmro2)/fs));
HbO = ctHb*(Sa*phi_a.*(1+va_t)+Sc*Far*phi_c+Sv*phi_v.*(1+vv_t) ...
    + (Sc/Sv*(Sc - Sv)*Far*phi_c*filter(b1,1,cbf - cmro2)/fs + (Sa-Sv)*phi_v*filter(b2,1,cbf - cmro2)/fs));
HbT = HbD + HbO;

dHbO = (HbO - O0);
dHbD = (HbD - D0);
dHbT = dHbO + dHbD;

Bold=(phi_a+phi_c+phi_v)*(3.4*(1-(HbD/D0))-((1-Sa)*va_t+(1-Sv)*vv_t)/(3-Sa-Sc-Sv));

figure, hold on;
% title('delta hemoglobin (uM)')
plot(T-90, dHbO, 'k--','LineWidth',2), plot(T-90, dHbD, 'k:','LineWidth',2), plot(T-90, dHbT, 'k-', 'LineWidth', 2);
legend('\DeltaHbO(t)','\DeltaHb(t)','\DeltaHbT(t)');
xlabel('Time (s)'), ylabel('Hemoglobin changes (\muM)')
set(gca, 'FontSize', 24, 'LineWidth', 2, 'Ylim', [-20 50]);
axis square
box off
xlim([20, 140])
set(gcf,'pos',[695   217   829   620]);
legend boxoff
savefig(gcf,'figures_new/hemodynamic_norm'), saveas(gcf,'figures_new/hemodynamic_norm.png')

figure, hold on;
plot(T-90, 100*HbO./HbT, 'k', 'LineWidth', 2);
xlabel('Time (s)'), ylabel('S_tO_2 (%)')
set(gca, 'FontSize', 24, 'LineWidth', 2, 'Ylim', [55 85]);
axis square
box off
xlim([20, 140])
set(gcf,'pos',[695   217   829   620]);
savefig(gcf,'figures_new/saturation_norm'), saveas(gcf,'figures_new/saturation_norm.png')

figure;
plot(T-90, Bold, 'k', 'LineWidth', 2);
xlabel('Time (s)'), ylabel('BOLD (a.u.)')
set(gca, 'FontSize', 24, 'LineWidth', 2, 'Ylim', [-0.04 0.08]);
axis square
box off
xlim([20, 140])
set(gcf,'pos',[695   217   829   620]);
savefig(gcf,'figures_new/BOLD_norm'), saveas(gcf,'figures_new/BOLD_norm.png')

figure, hold on;
plot(T-90, cbf,'k','LineWidth', 2), plot(T-90, cbv, 'k--','LineWidth', 2), plot(T-90, cmro2, 'k:','LineWidth', 2);
legend('cbf', 'cbv', 'cmro2')
xlabel('Time (s)'), ylabel('Fractional change (a.u.)')
set(gca, 'FontSize', 24, 'LineWidth', 2, 'Ylim', [-0.5 1.5]);
axis square
box off
xlim([20, 140])
set(gcf,'pos',[695   217   829   620]);
legend boxoff
savefig(gcf,'figures_new/cbf_norm'), saveas(gcf,'figures_new/cbf_norm.png')

