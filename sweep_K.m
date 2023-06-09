clear;

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

t = 0:0.01:60;
INIT = [n, Vm, Ca, 0., 0.5, 0.7, 0.1709];
% Enumerate through a number of Kp and CPP levels
conc_K = 2:0.2:40;
CPP = 50:10:100;
R = zeros(length(t), length(conc_K), length(CPP));

Opt = odeset('MaxStep', 0.05, 'BDF', 'on');
for j = 1:length(CPP)
    for i=1:length(conc_K)
        [~, result] = ode15s(@(t,y)vasodynamic(t, y, conc_K(i), CPP(j)), t, INIT);
        R(:,i,j) = result(:,end)/2/pi;
    end
end

avgR = squeeze(mean(R(1001:end, :, :), 1));
stdR = squeeze(std(R(1001:end, :, :), [], 1));
% save('sweep_K')

%% Generating Fig. 8(a)
% Set line colors for each CPP level
linecolors = [0.94,0.94,0.94;
0.83,0.82,0.78;
0.80,0.80,0.80;
0.50,0.50,0.50;
0.31,0.31,0.31;
0.00,0.00,0.00];

figure, hold on;
str_legend = split(num2str(CPP));
set(gcf,'pos',[140,32,1142,822])
h_patch = [];
for i=1:length(CPP)
    h=shadedErrorBar(conc_K, avgR(:,i), stdR(:,i),'lineProps',{'color',linecolors(i,:)},'patchSaturation',0.5);
    h_patch = [h_patch; h.patch];
    str_legend{i} = [str_legend{i}, 'mmHg'];
end
axis square
xlim([2,40])
xlabel('Perivascular K^+ Concentration (mM)'), ylabel({'Average vessel radius', 'with standard deviation (cm)'})
legend(flipud(h_patch), flipud(str_legend))
legend boxoff
set(gca, 'FontSize', 40);

%% Generating Fig. 8(b)
figure, hold on
set(gcf,'pos',[125,28,1051,866])
plot(t, squeeze(R(:,5,2)), 'k-', 'LineWidth',3)
plot(t, squeeze(R(:,5,6)), 'k:', 'LineWidth',3)
xlim([10,60]),ylim([0.021,0.026]), xlabel('Time (s)'), ylabel('Radius (cm)')
axis square
set(gca, 'FontSize', 40)
legend('60 mmHg', '100 mmHg','Location','east')
