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
str_legend = split(num2str(CPP));
figure, hold on;
for i=1:length(CPP)
    errorbar(conc_K, avgR(:,i), stdR(:,i));
    str_legend{i} = [str_legend{i}, 'mmHg'];
end
set(gca, 'FontSize', 18);
xlabel('Perivascular K^+ Concentration'), ylabel('Mean Vessel Radius with Standard Error (cm)')
legend(str_legend)

save('sweep_K')
