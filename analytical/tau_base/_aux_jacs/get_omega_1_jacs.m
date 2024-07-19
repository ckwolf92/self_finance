% type 1 matrices

omega = omega_1;

% auxiliary matrix

A_1 = zeros(2*T,2*T); % order (c,d) & (BC, EE)
for t = 1:T
    A_1(t,t) = 1;
    A_1(t,T+t) = 1;
    if t > 1
        A_1(t,T+t-1) = -1/beta;
    end
end
for t = T+1:2*T
    if t == T+1
        A_1(t,t-T) = 1 - omega_1 * (1-beta*omega_1);
        A_1(t,t-T+1) = - beta * omega_1;
    elseif t < 2*T
        A_1(t,t-T) = 1 - omega_1 * (1-beta*omega_1);
        A_1(t,t-T+1) = - beta * omega_1;
        A_1(t,t-1) = - (1-beta*omega_1) * (1-omega_1) * 1/beta;
    elseif t == 2*T
        A_1(t,t-T) = 1 - omega_1 * (1-beta*omega_1);
        A_1(t,t-1) = 1;
    end
end

A_1_inv = A_1^(-1);

% baseline

[c_base,d_base] = c_olg_fn(exo,A_1_inv);

% income

C_y_1 = NaN(T,T);
D_y_1 = NaN(T,T);

for t = 1:T
    exo.y_hat(t) = step;
    [c_shock,d_shock] = c_olg_fn(exo,A_1_inv);
    C_y_1(:,t) = (c_shock-c_base)/step;
    D_y_1(:,t) = (d_shock-d_base)/step;
    exo.y_hat(t) = 0;
end

% nominal interest rates

C_i_1 = NaN(T,T);
D_i_1 = NaN(T,T);

for t = 1:T
    exo.i_hat(t) = step;
    [c_shock,d_shock] = c_olg_fn(exo,A_1_inv);
    C_i_1(:,t) = (c_shock-c_base)/step;
    D_i_1(:,t) = (d_shock-d_base)/step;
    exo.i_hat(t) = 0;
end

% inflation

C_pi_1 = NaN(T,T);
D_pi_1 = NaN(T,T);

for t = 1:T
    exo.pi_hat(t) = step;
    [c_shock,d_shock] = c_olg_fn(exo,A_1_inv);
    C_pi_1(:,t) = (c_shock-c_base)/step;
    D_pi_1(:,t) = (d_shock-d_base)/step;
    exo.pi_hat(t) = 0;
end

% demand shock

C_d_1 = NaN(T,T);
D_d_1 = NaN(T,T);

for t = 1:T
    exo.zeta_hat(t) = step;
    [c_shock,d_shock] = c_olg_fn(exo,A_1_inv);
    C_d_1(:,t) = (c_shock-c_base)/step;
    D_d_1(:,t) = (d_shock-d_base)/step;
    exo.zeta_hat(t) = 0;
end

omega = [];