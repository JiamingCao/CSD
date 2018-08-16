function dydt = vasodynamic(t, y, Kp, P)
    n = y(1);  % NO UNIT
    Vm = y(2);  % mV
    Ca_smc = y(3);  % nM
    k = y(4);  % NO UNIT
    w = y(5);   % NO UNIT
    yprime = y(6);  % NO UNIT
    x = y(7);   % cm
     %% Constants
    g_leak = 7.85e-14;  % C s-1 mV-1 (ie g_L in SMC model)
    v_leak = -70.;  % mV (ie v_L in SMC model)
    v_ca = 80.0;  % mV
    v_K = -90.0;  % mV
    v1 = -17.4 - 0.06 * P;    % mV, P in mmHg
    v2 = 25.0;  % mV
    v4 = 14.5;  % mV
    v5 = 8.0;  % mV
    v6 = -15.0;  % mV
    Ca3 = 400.;  % nM
    Ca4 = 150.;  % nM
    phi_n = 2.664;  % s-1

    g_ca = 1.57e-13;  % C s-1 mV-1
    g_K_smc = 3.14e-13;  % C s-1 mV-1
    Kd = 1e3;  % nM
    BT = 1e5;  % nM
    alpha = 1.1996e16;  % nMC-1,originally 8e-15
    k_ca = 1.36e2;  % s-1
    C_smc = 1.4687e-14;  % C mV-1
    %% Mechanics constants
    x1 = 1.2;
    x2 = 1.3e-1;
    x3 = 2.22;
    x4 = 7.11e-1;
    x5 = 8.0e-1;
    x6 = 1.0e-2;
    x7 = 3.21e-1;
    x8 = 8.9e-1;
    x9 = 9.05e-3;
    u1 = 41.76;
    u2 = 4.73e-2;
    u3 = 5.84e-2;
    y0 = 9.28e-1;
    y1 = 6.39e-1;
    y2 = 3.5e-1;
    y3 = 7.88e-1;
    y4 = 8.0e-1;
    q = 3;
    psi_m = 0.3;
    a = 0.281;
    b = 5.0;
    c = 3e-2;
    d = 1.3;
    Ca_m = 3.5e2;  % nM
    Ca_ref = 4e2;  % nM
    sigma_y0_sharp = 2.6e6;  % dyne cm-2
    sigma0_sharp = 3e6;  % dyne cm-2
    k_psi = 3.3;  % s-1
    v_ref = 0.24;  % s-1
    x0 = 0.213;  % cm
    S = 9.0e-3;  % cm2
    A = 1.5e-3;  % cm2
    w_e = 0.9;
    w_m = 0.7;
    tau = 0.2;  % dyne cm-1 s
    delta_p = P * 1.3328e3;     % dyne cm-2
    %% The model
    dydt = zeros(7,1);
    
    v3 = -0.5 * v5 * tanh((Ca_smc - Ca3) / Ca4) + v6;
    lambda_n = phi_n * cosh((Vm - v3) / 2 / v4);
    n_infty = 0.5 * (1 + tanh((Vm - v3) / v4));
    dydt(1) = lambda_n * (n_infty - n);

    % KIR Enters here
    v_KIR = 57 * log10(Kp) - 130;  % mV
    g_KIR = 145 * sqrt(Kp) * 1e-15;  % C s-1 mV-1
    % End KIR

    m_infty = 0.5 * (1 + tanh((Vm - v1) / v2));
    dydt(2) = -(g_leak * (Vm - v_leak) + g_K_smc * n * (Vm - v_K) + g_ca * m_infty * (Vm - v_ca) + g_KIR * k * (Vm - v_KIR)) / C_smc;

    rho = power(Kd + Ca_smc, 2) / (power(Kd + Ca_smc, 2) + Kd * BT);
    dydt(3) = -(alpha * g_ca * m_infty * (Vm - v_ca) + k_ca * Ca_smc) * rho;

    alpha_KIR = 1020. / (1. + exp((Vm - v_KIR + 18.) / 5.5)); % Farr&David, 6.8
%     beta_KIR = 26.9 * exp(0.06 * (Vm - v_KIR + 18.));
    beta_KIR = 26.9 * exp(0.06 * (Vm - v_KIR + 29));    % Farr&David says 18, Kurashi originally says 29
    dydt(4) = (alpha_KIR + beta_KIR) * (alpha_KIR / (alpha_KIR + beta_KIR) - k);

    psi = func_psi(Ca_smc, Ca_m, q);
    dydt(5) = k_psi * (psi / (psi + psi_m) - w);

    u_prime = x / x0 - yprime;
    sigma_u = u2 * exp(u1 * u_prime) - u3;
    w_ref = func_psi(Ca_ref, Ca_m, q) / (psi_m + func_psi(Ca_ref, Ca_m, q));
    sigma_y0 = sigma_y0_sharp / w_ref * w;
    sigma_y = (sigma_y0 / sigma0_sharp) * (exp(-0.5 * power(yprime - y0, 2) / power(y1 / (yprime + y2), 2. * y4)) - y3) / (1. - y3);
    psi_ref = func_psi(Ca_ref, Ca_m, q);
    tmp = sigma_u / sigma_y;
    if tmp < 1.
        dydt(6) = -v_ref * (psi / psi_ref) * a * (1. - tmp) / (a + tmp);
    else
        dydt(6) = c * (exp(b * (tmp - d)) - exp(b * (1. - d)));
    end

    if x > sqrt(pi * A)
        f_delta_p = 0.5 * delta_p * (x / pi - A / x);
    else
        f_delta_p = 0.;
    end

    x_prime = x / x0;
    sigma_x = x3 * (1 + tanh((x_prime - x1) / x2)) + x4 * (x_prime - x5) - x8 * power(x6 / (x_prime - x7), 2) - x9;
    f_x = w_e * S * sigma_x * sigma0_sharp;
    f_u = w_m * S * sigma_u * sigma0_sharp;

    dydt(7) = (f_delta_p - f_x - f_u) / tau;
     function psi = func_psi(ca, ca_m, q)
        psi = ca^q/(ca^q + ca_m^q);