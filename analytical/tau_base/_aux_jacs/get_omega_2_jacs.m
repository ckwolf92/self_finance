% type 2 matrices

omega = omega_2;

% auxiliary matrix

A_2 = zeros(2*T,2*T); % order (c,d) & (BC, EE)
for t = 1:T
    A_2(t,t) = 1;
    A_2(t,T+t) = 1;
    if t > 1
        A_2(t,T+t-1) = -1/beta;
    end
end
for t = T+1:2*T
    if t == T+1
        A_2(t,t-T) = 1 - omega_2 * (1-beta*omega_2);
        A_2(t,t-T+1) = - beta * omega_2;
    elseif t < 2*T
        A_2(t,t-T) = 1 - omega_2 * (1-beta*omega_2);
        A_2(t,t-T+1) = - beta * omega_2;
        A_2(t,t-1) = - (1-beta*omega_2) * (1-omega_2) * 1/beta;
    elseif t == 2*T
        A_2(t,t-T) = 1 - omega_2 * (1-beta*omega_2);
        A_2(t,t-1) = 1;
    end
end

A_2_inv = A_2^(-1);

% baseline

[c_base,d_base] = c_olg_fn(exo,A_2_inv);

% income

C_y_2 = NaN(T,T);
D_y_2 = NaN(T,T);

for t = 1:T
    exo.y_hat(t) = step;
    [c_shock,d_shock] = c_olg_fn(exo,A_2_inv);
    C_y_2(:,t) = (c_shock-c_base)/step;
    D_y_2(:,t) = (d_shock-d_base)/step;
    exo.y_hat(t) = 0;
end

% nominal interest rates

C_i_2 = NaN(T,T);
D_i_2 = NaN(T,T);

for t = 1:T
    exo.i_hat(t) = step;
    [c_shock,d_shock] = c_olg_fn(exo,A_2_inv);
    C_i_2(:,t) = (c_shock-c_base)/step;
    D_i_2(:,t) = (d_shock-d_base)/step;
    exo.i_hat(t) = 0;
end

% inflation

C_pi_2 = NaN(T,T);
D_pi_2 = NaN(T,T);

for t = 1:T
    exo.pi_hat(t) = step;
    [c_shock,d_shock] = c_olg_fn(exo,A_2_inv);
    C_pi_2(:,t) = (c_shock-c_base)/step;
    D_pi_2(:,t) = (d_shock-d_base)/step;
    exo.pi_hat(t) = 0;
end

% demand shock

C_d_2 = NaN(T,T);
D_d_2 = NaN(T,T);

for t = 1:T
    exo.zeta_hat(t) = step;
    [c_shock,d_shock] = c_olg_fn(exo,A_2_inv);
    C_d_2(:,t) = (c_shock-c_base)/step;
    D_d_2(:,t) = (d_shock-d_base)/step;
    exo.zeta_hat(t) = 0;
end

omega = [];