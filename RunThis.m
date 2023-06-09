clear;
%% Tuckwell model
percent = 100;  % Reduce efficiency of K and neurotransmitter cleaning up; 
% 100, 85, or 35 for mild, intermediate and severe CSD (see Results)
[T_vec, K_result, Ca_result, Na_result, Cl_result, glu_result, gaba_result] = genTuckwell2D(percent/100);
save(['tuckwell_',num2str(percent),'percent'], 'Ca_result', 'Cl_result', 'K_result', 'Na_result', 'gaba_result', 'glu_result', 'T_vec','percent');

% Generating the two subplots in Fig. 3. This is done when percent=100
figure, mesh(linspace(-5.2, 5.2, 201), linspace(-5.2, 5.2, 201), K_result(:,:,463))
axis tight,colormap gray,view(-45,75),colorbar
xlabel('x (mm)'), ylabel('y (mm)'), zlabel('Ks (mM)')
set(gca,'FontSize',18)

figure, mesh(linspace(-5.2, 5.2, 201), linspace(-5.2, 5.2, 201), glu_result(:,:,463))
axis tight,colormap gray,view(-45,75),colorbar
xlabel('x (mm)'), ylabel('y (mm)'), zlabel('Glu (mM)')
set(gca,'FontSize',18)
%% Neurovascular coupling model
% First add some baseline before KCl application
dt = T_vec(2) - T_vec(1);
T = 26*[0:dt:3.5-dt, 3.5+T_vec]'; % Scale time axis. One unit time is 26 sec
K = cat(3, K_result(1)*ones([201,201,3.5/dt]), K_result);
Glu = cat(3, glu_result(1)*ones([201,201,3.5/dt]), glu_result);

% Probe at 1.56 mm away from the center
K_probe = squeeze(K(101, 131, :));
Glu_probe = squeeze(Glu(101,131,:));

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
INIT = [Ca_astr, Ca_ER, IP3, h, EET, n_BK, 100., Vk, n, Vm, Ca, 0., 5., 0.5, 0.7, 0.1709];% set up initial values
% Now solve the NVC model using built-in PDE solver
serca_ratio = 1;    % Efficiency of SERCA pump from 0 to 1 (0% - 100%)
% 1 for mild, 0.9 for intermediate, and 0.5 for severe
BK_shift = 0;    % Negative shift of v6,BK to adjust the sensitity of BK channels
% 0 for mild, -1 for intermediate, and -12 for severe
opt=odeset('MaxStep',0.1,'BDF','on');
[T,result]=ode15s(@(t,y)nvcoupling(t, y, K_probe, Glu_probe, T, 80, serca_ratio, BK_shift),t,INIT,opt); % CPP=80mmHg; the last term is 1 or 0.5 for full or 50% SERCA pump power
R = result(:,end)/2/pi;
Ca = result(:,1);
Kp = result(:, 13);
save(['ProbedRadius_Uptake', num2str(percent), '_SERCA', num2str(100*serca_ratio), '_BKshift', num2str(BK_shift)], 'T', 'R', 'Kp','Ca');

%% Generating the subplots in Fig. 4. 
% (a) When percent=100, serca_ratio=1;
% (b) when percent=100, serca_ratio=0.5
figure;
pos = [ 471, 97, 1105, 867];
set(gcf, 'pos', pos, 'defaultAxesColorOrder', [0,0,0; 0,0,0]);
plot(T-90, Ca, 'k--', 'LineWidth', 2), ylabel('Ca_{astr}(nM)');
set(gca, 'FontSize', 32);
ylim([120 600])
yyaxis right
plot(T-90, Kp, 'k-', 'LineWidth', 2), ylabel('Kp (mM)')
set(gca, 'FontSize', 32);
xlim([0 150])
ylim([3 60])
box off
legend('Ca','Kp')
legend boxoff
xlabel('Time (s)')

%% Generating subplots in Fig. 5. 
% (a) when percent=100,serca_ratio=1,BK_shift=0
% (b) when percent=35 and the other two the same
T2 = 26*[0:dt:3.5-dt, 3.5+T_vec]';
figure, set(gcf, 'pos', pos), hold on;
plot(T-90, interp1(T2, K_probe, T), 'k-', 'LineWidth', 2);
plot(T-90, Kp, 'k--', 'LineWidth', 2);
plot(T-90, interp1(T2, Glu_probe, T), 'k:', 'LineWidth', 2);
xlim([0,150]), ylim([0, 55]);
set(gca, 'FontSize', 32, 'LineWidth', 2);
xlabel('Time (s)'), ylabel('Concentration (mM)')
legend('Ks', 'Kp', 'Glu', 'Ca'), legend boxoff

%% Generating subplots in Fig. 6. 
% (a) Everything normal value 
% (b) only BK_shift=-12, everything else the same
figure, set(gcf, 'pos', pos), hold on;
plot(T-90, interp1(T2, K_probe, T), 'k-', 'LineWidth', 2);
plot(T-90, Kp, 'k--', 'LineWidth', 2);
xlim([0,150]), ylim([0, 70]);
set(gca, 'FontSize', 32, 'LineWidth', 2);
xlabel('Time (s)'), ylabel('Concentration (mM)')
legend('Ks', 'Kp'), legend boxoff

%% Generating subplots in Fig. 7
% (a) BK_shift=0, other parameters same as mild case
% (b) BK_shift=-1, other parameters same as mild case
% (c) BK_shift=-12, other parameters same as mild case
pos = [293    60   892   806];
figure, set(gcf, 'pos', pos, 'defaultAxesColorOrder', [0,0,0; 0,0,0; 0,0,0]);
subplot(2,1,1), plot(T-90, Kp, 'k', 'LineWidth', 2);
ax = gca;
ax.XAxis.Visible = 'off';
ylabel('Kp (mM)')
xlim([0, 150]), ylim([0,70]);
yticks(0:10:70);
set(ax, 'Fontsize', 24);
box off

subplot(2,1,2), plot(T-90, R, 'k', 'LineWidth', 2);
ylabel('Vessel Radius (cm)'), xlabel('Time (s)')
xlim([0, 150]), ylim([0.015,0.03]);
yticks(0.015:0.005:0.03);
set(gca, 'Fontsize', 24);
box off
yyaxis right
R_baseline = mean(R(500:1500));
plot(T-90,100*(R - R_baseline)/R_baseline,'Color','None')
ylabel('% Change'), ylim([100*(0.015 - R_baseline)/R_baseline, 100*(0.03 - R_baseline)/R_baseline])
% yticks(-20:10:20)   % Use this for subplots (a) (b)
yticks(-40:20:0)    % Use this for subplot (c)
box off
set(gca, 'Fontsize', 24);

%% Run the hemodynamics model
T0=100; %%% uM
Tc      = 0.3;  %%% capillary transit time
Tv      = 6;    %%% venous transit time
pa      = 0.04;        %%% CBV_a_0 / CBV_0 ... arterial contribution. 0.04 means 4%
pc      = 0.55;  %%% F*CBV_c_0 / CBV_0 ... capillary contribution with Fahreaus factor 
alpha   = 1.2;
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

cbf = (cbv + 1).^(2.6)-1;

% Get metabolism using balloon model
E0 = 0.4;
cmro2 = (1-(1-E0).^(1./(1+cbf)))/E0.*(1+cbf) - 1;

% Fantini model (CHS paper)
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

%% Generating subplots in Fig. 9, Fig. 10, and Supplementary Fig. 1
% Fig. 9 (a)-(c): percent=100, serca_ratio=1, BK_shift=0
% Fig. 9 (d)-(f): percent=100, serca_ratio=0.5, BK_shift=0
% Fig. 9 (g)-(i): percent=35, serca_ratio=1, BK_shift=0
% Fig. 9 (j)-(l): percent=100, serca_ratio=1, BK_shift=-12
% Fig. 10 (a), Supp. Fig. 1 (a)(d): percent=100, serca_ratio=0.9, BK_shift=0
% Fig. 10 (b), Supp. Fig. 1 (b)(e): percent=85, serca_ratio=1, BK_shift=0
% Fig. 10 (c), Supp. Fig. 1 (c)(f): percent=100, serca_ratio=1, BK_shift=-1

figure, hold on;
plot(T-90, dHbO, 'k--','LineWidth',2), plot(T-90, dHbD, 'k:','LineWidth',2), plot(T-90, dHbT, 'k-', 'LineWidth', 2);
legend('\DeltaHbO(t)','\DeltaHb(t)','\DeltaHbT(t)');
xlabel('Time (s)'), ylabel('Hemoglobin changes (\muM)')
set(gca, 'FontSize', 24, 'LineWidth', 2, 'Ylim', [-25 35]);
axis square
box off
xlim([20, 140])
set(gcf,'pos',[695   217   829   620]);
legend boxoff
% savefig(gcf,'figures_new2/hemodynamic_mid_Ca'), saveas(gcf,'figures_new2/hemodynamic_mid_Ca.png')

figure;
plot(T-90, Bold*100, 'k', 'LineWidth', 2);
xlabel('Time (s)'), ylabel('BOLD (%)')
set(gca, 'FontSize', 24, 'LineWidth', 2, 'Ylim', [-2 5]);
axis square
box off
xlim([20, 140])
set(gcf,'pos',[695   217   829   620]);
% savefig(gcf,'figures_new2/BOLD_mid_Ca'), saveas(gcf,'figures_new2/BOLD_mid_Ca.png')

figure, hold on;
plot(T-90, cbf,'k','LineWidth', 2), plot(T-90, cbv, 'k--','LineWidth', 2), plot(T-90, cmro2, 'k:','LineWidth', 2);
legend('cbf', 'cbv', 'cmro2')
xlabel('Time (s)'), ylabel('Fractional change (a.u.)')
set(gca, 'FontSize', 24, 'LineWidth', 2, 'Ylim', [-0.65 1.5]);
axis square
box off
xlim([20, 140])
set(gcf,'pos',[695   217   829   620]);
legend boxoff
% savefig(gcf,'figures_new2/cbf_mid_Ca'), saveas(gcf,'figures_new2/cbf_mid_Ca.png')