function dydt = nvcoupling(t, y, inputK, inputGlu, inputT, P, ratioPump)
 % Modified David astrocyte assuming that potential change comes only from flux of potassium
    Ca_astr = y(1);  % nM
    Ca_ER = y(2);  % nM
    IP3 = y(3);  % uM
    h_k = y(4);  % NO UNIT
    EET = y(5);  % uM
    n_BK = y(6);  % NO UNIT
    K_i = y(7);  % mM
    Vk = y(8);  % mV
    n = y(9);  % NO UNIT
    Vm = y(10);  % mV
    Ca_smc = y(11);  % nM
    k = y(12);  % NO UNIT
    Kp = y(13);  % mM
    w = y(14);   % NO UNIT
    yprime = y(15);  % NO UNIT
    x = y(16);   % cm

    %% Constants
    J_max = 2.88e6;  % nM s-1
    K_I = 0.03;  % uM
    K_act = 170;  % nM
    V_max = 2e4;  % nM s-1
    k_pump = 240;  % uM
    P_L = 5.2e3;  % nM s-1 (Fink says 0.0804 uM s-1)
    B_cyt = 2.44e-2;  % NO UNIT
    VR_ERcyt = 0.185;  % NO UNIT
    delta = 1.235e-3;  % NO UNIT
    K_G = 8.82;  % NO UNIT
    % r_h = 0.8           % uM s-1
    r_h = 5.2;   % uM s-1
    k_deg = 1.25;  % s-1
    % k_on = 2.0          % uM-1 s-1 (Fink says 8.0)
    k_on = 8.0;  % uM-1 s-1
    K_inh = 0.1;  % uM
    V_eet = 72.;  % s-1
    k_eet = 7.2;  % s-1
    Ca_min = 0.1;  % uM
    g_leak = 7.85e-14;  % C s-1 mV-1 (ie g_L in SMC model)
    g_leak_astr = 6.67e-12; % C s-1 mV-1
    v_leak = -70.;  % mV (ie v_L in SMC model)
    g_BK = 1.5e-11;  % C s-1 mV-1
    v_BK = -95.;  % mV
    v_ca = 80.0;  % mV
    v_K = -90.0;  % mV
%     J_NaKmax_prime = 15e-11 * ratioPump;  % C s-1
%     J_NaKmax_prime = 1e-10 * ratioPump;  % C s-1
    J_NaKmax_prime = 1.8e-10 * ratioPump;  % C s-1
    KKoa = 1.5;  % mM
    C_astr = 4e-14;  % C mV-1
    EET_shift = 2.;  % mV uM-1
%     v1 = -22.5;  % mV
    v1 = -17.4 - 0.06 * P;    % mV, P in mmHg
    v2 = 25.0;  % mV
    v4 = 14.5;  % mV
    v5 = 8.0;  % mV
    v6 = -15.0;  % mV
    Ca3 = 400.;  % nM
    Ca4 = 150.;  % nM
    phi_n = 2.664;  % s-1

    g_ca = 1.57e-13;  % C s-1 mV-1
    g_k_astr = 2.82e-11;  % C s-1 mV-1
%     g_k_astr = 5.67e-11;  % C s-1 mV-1
    g_K_smc = 3.14e-13;  % C s-1 mV-1
    g_nkcc = 4e-14;  % C s-1 mV-1
    g_kcc = 1.4e-12;  % C s-1 mV-1
    F = 9.649e4;  % C mol-1
    R = 8.31;  % V C mol-1 K-1
    T = 310;  % K
    Na_i = 15;  % mM
    Cl_i = 25;  % mM
    Cl_o = 145;  % mM
    Na_o = 146;  % mM
    Kd = 1e3;  % nM
    BT = 1e5;  % nM
    alpha = 1.1996e16;  % nMC-1,originally 8e-15
    k_ca = 1.36e2;  % s-1
    Kp_min = 3.;  % mM
    C_smc = 1.4687e-14;  % C mV-1
%     C_smc = 1.96e-14;  % C mV-1, original

    %% Perivascular constants
    VR_pa = 8e-4;  % NO UNIT
    VR_ps = 0.1;  % NO UNIT
    R_decay = 0.1;  % s-1

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
    % delta_p = 8e4           % dyne cm-2
    delta_p = P * 1.3328e3;     % dyne cm-2
    

    %% Dummy input
    % Ks = dummy_Ks(t)
    % if Ks < 4.
    %     rho_G = 0.
    % else
    %     rho_G = Ks / 60.

    %% Better input
    Ks = get_K(inputK, inputT, t);
%     Cl_o = get_K(inputCl, inputT, t);
%     Na_o = get_K(inputNa, inputT, t);
    rho_G = get_rho(inputGlu, inputT, t);

    %% The model
    dydt = zeros(16,1);

    J_NaK = J_NaKmax_prime * Ks / (Ks + KKoa);
    J_IP3 = J_max * ((IP3 / (IP3 + K_I) * Ca_astr / (Ca_astr + K_act) * h_k) ^ 3) * (1 - Ca_astr / Ca_ER);
    J_pump = V_max * (Ca_astr ^ 2 / (k_pump ^ 2 + Ca_astr ^ 2));
    J_ERleak = P_L * (1 - Ca_astr / Ca_ER);
    dydt(1) = B_cyt * (J_IP3 - J_pump + J_ERleak);

    dydt(2) = -(1 / VR_ERcyt) * dydt(1);

    G = (rho_G + delta) / (K_G + rho_G + delta);
    dydt(3) = r_h * G - k_deg * IP3;

    dydt(4) = k_on * (K_inh - (Ca_astr / 1000 + K_inh) * h_k);

    dydt(5) = V_eet * (Ca_astr / 1000 - Ca_min) - k_eet * EET;

    v3_astr = -0.5 * v5 * tanh((Ca_astr - Ca3) / Ca4) + v6;
    lambda_BK = phi_n * cosh((Vk - v3_astr) / 2 / v4);
    n_BK_infty = 0.5 * (1 + tanh((Vk + EET_shift * EET - v3_astr) / v4));
    dydt(6) = lambda_BK * (n_BK_infty - n_BK);

%     J_NKCC = 1000 * g_nkcc * log(Ks * Cl_o / K_i / Cl_i) * R * T / F;
    J_KCC = 1000 * g_kcc * log(Ks * Cl_o / K_i / Cl_i) * R * T / F;
%     J_KCC = 1000 * g_kcc * log(Na_o * Ks * Cl_o ^ 2 / K_i / Na_i / Cl_i ^ 2) * R * T / F;
    J_NKCC = 1000 * g_nkcc * log(Na_o * Ks * Cl_o ^ 2 / K_i / Na_i / Cl_i ^ 2) * R * T / F;
    dydt(7) = (-g_k_astr * (Vk - (v_K)) + 2 * J_NaK - g_BK * n_BK * (Vk - v_BK) + J_NKCC + J_KCC) / 5.79e-11;
%     dydt(7) = ( -g_k_astr * (Vk - (v_K)) + 2*J_NaK - g_BK * n_BK * (Vk - v_BK) + J_NKCC + J_KCC) / 5.79e-11;
    dydt(8) = (-g_k_astr * (Vk - (v_K)) + 2 * J_NaK - g_BK * n_BK * (Vk - v_BK) + J_NKCC + J_KCC) / C_astr;
%     dydt(8) = (-g_leak*(Vk-v_leak) -g_k_astr * (Vk - (v_K)) + J_NaK - g_BK * n_BK * (Vk - v_BK) + J_NKCC + J_KCC) / C_astr;

    v3 = -0.5 * v5 * tanh((Ca_smc - Ca3) / Ca4) + v6;
    lambda_n = phi_n * cosh((Vm - v3) / 2 / v4);
    n_infty = 0.5 * (1 + tanh((Vm - v3) / v4));
    dydt(9) = lambda_n * (n_infty - n);

    % KIR Enters here
    v_KIR = 57 * log10(Kp) - 130;  % mV
    g_KIR = 145 * sqrt(Kp) * 1e-15;  % C s-1 mV-1
    % End KIR

    m_infty = 0.5 * (1 + tanh((Vm - v1) / v2));
    dydt(10) = -(g_leak * (Vm - v_leak) + g_K_smc * n * (Vm - v_K) + g_ca * m_infty * (Vm - v_ca) + g_KIR * k * (Vm - v_KIR)) / C_smc;

    rho = power(Kd + Ca_smc, 2) / (power(Kd + Ca_smc, 2) + Kd * BT);
    dydt(11) = -(alpha * g_ca * m_infty * (Vm - v_ca) + k_ca * Ca_smc) * rho;

    alpha_KIR = 1020. / (1. + exp((Vm - v_KIR + 18.) / 5.5)); % Farr&David, 6.8, 5.5 was used
%     beta_KIR = 26.9 * exp(0.06 * (Vm - v_KIR + 18.));
    beta_KIR = 26.9 * exp(0.06 * (Vm - v_KIR + 29));    % Farr&David says 18, Kurashi originally says 29
    dydt(12) = (alpha_KIR + beta_KIR) * (alpha_KIR / (alpha_KIR + beta_KIR) - k);

    I_BK = g_BK * n_BK * (Vk - v_BK);
    J_BK = I_BK / (5.79e-11);
    J_KIR = g_KIR * k * (Vm - v_KIR) / (5.79e-11);
    dydt(13) = J_BK / VR_pa + J_KIR / VR_ps - R_decay * (Kp - Kp_min);

    psi = func_psi(Ca_smc, Ca_m, q);
    dydt(14) = k_psi * (psi / (psi + psi_m) - w);

    u_prime = x / x0 - yprime;
    sigma_u = u2 * exp(u1 * u_prime) - u3;
    w_ref = func_psi(Ca_ref, Ca_m, q) / (psi_m + func_psi(Ca_ref, Ca_m, q));
    sigma_y0 = sigma_y0_sharp / w_ref * w;
    sigma_y = (sigma_y0 / sigma0_sharp) * (exp(-0.5 * power(yprime - y0, 2) / power(y1 / (yprime + y2), 2. * y4)) - y3) / (1. - y3);
    psi_ref = func_psi(Ca_ref, Ca_m, q);
    tmp = sigma_u / sigma_y;
    if tmp < 1.
        dydt(15) = -v_ref * (psi / psi_ref) * a * (1. - tmp) / (a + tmp);
    else
        dydt(15) = c * (exp(b * (tmp - d)) - exp(b * (1. - d)));
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

    dydt(16) = (f_delta_p - f_x - f_u) / tau;
    
    function K = get_K(K,T,t)
        K = interp1(T,K,t,'linear');

    function rho = get_rho(glu, T, t)
        g = interp1(T,glu,t,'linear');
        rho=g/(g+2);
    
    function psi = func_psi(ca, ca_m, q)
        psi = ca^q/(ca^q + ca_m^q);