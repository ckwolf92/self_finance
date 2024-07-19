% inflation

pi_seq = zeros(T,1);
pi_seq(T) = kappa * y_seq(T) - kappa/2 * g_seq(T);
for t = T-1:-1:1
    pi_seq(t) = kappa * y_seq(t) + beta * pi_seq(t+1) - kappa/2 * g_seq(t); % \sigma = \varphi = 1
end

% monetary block

i_seq = [pi_seq(2:T);0] + phi * y_seq;

% transfers and debt

d_seq     = NaN(T,1);
t_seq     = NaN(T,1);

% H rules

if isempty(tau_d) == 1 % this will be the case for the H rules

if H == 0

d_seq(1) = 0;
t_seq(1) = g_seq(1) + 1/beta * 0 + D_SS/Y_SS * (0 - pi_seq(1)) - d_seq(1);

for t = 2:T
    d_seq(t) = 0;
    t_seq(t) = g_seq(t) + 1/beta * d_seq(t-1) + D_SS/Y_SS * (i_seq(t-1) - pi_seq(t)) - d_seq(t);
end

else

t_seq(1) = tau_y * y_seq(1);
d_seq(1) = 1/beta * 0 - t_seq(1) + g_seq(1) + D_SS/Y_SS * (0 - pi_seq(1));

if H > 1

for t = 2:H
    t_seq(t)     = tau_y * y_seq(t);
    d_seq(t)     = 1/beta * d_seq(t-1) - t_seq(t) + g_seq(t) + D_SS/Y_SS * (i_seq(t-1)-pi_seq(t));
end

end

for t = H+1:T
    d_seq(t) = 0;
    t_seq(t) = g_seq(t) + 1/beta * d_seq(t-1) + D_SS/Y_SS * (i_seq(t-1) - pi_seq(t)) - d_seq(t);
end

end

end

% tau_d rules

if isempty(tau_d) == 0
    
if H == 0

d_seq(1) = 0;
t_seq(1) = g_seq(1) + 1/beta * 0 + D_SS/Y_SS * (0 - pi_seq(1)) - d_seq(1);

for t = 2:T
    d_seq(t) = 0;
    t_seq(t) = g_seq(t) + 1/beta * d_seq(t-1) + D_SS/Y_SS * (i_seq(t-1) - pi_seq(t)) - d_seq(t);
end

else

t_seq(1) = tau_d * (1/beta * 0 + D_SS/Y_SS * (0 - pi_seq(1))) + tau_y * y_seq(1) + tau_d * (1-tau_y) * g_seq(1);
d_seq(1) = 1/beta * 0 - t_seq(1) + g_seq(1) + D_SS/Y_SS * (0 - pi_seq(1));

if H > 1

for t = 2:H
    t_seq(t)     = tau_y * y_seq(t) + tau_d * (1-tau_y) * g_seq(t);
    d_seq(t)     = 1/beta * d_seq(t-1) - t_seq(t) + g_seq(t) + D_SS/Y_SS * (i_seq(t-1)-pi_seq(t));
end

end

for t = H+1:T
    t_seq(t)     = tau_d * ((1-tau_y) * g_seq(t) + 1/beta * d_seq(t-1) + D_SS/Y_SS * (i_seq(t-1) - pi_seq(t))) + tau_y * y_seq(t);
    d_seq(t)     = 1/beta * d_seq(t-1) - t_seq(t) + g_seq(t) + D_SS/Y_SS * (i_seq(t-1) - pi_seq(t));
end

end

end

% consumption

c_seq    = C_y * (y_seq - t_seq) + C_i * i_seq + C_pi * pi_seq;
d_seq_HH = D_y * (y_seq - t_seq) + D_i * i_seq + D_pi * pi_seq;