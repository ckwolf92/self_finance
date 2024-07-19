% inflation

pi_seq = zeros(T,1);
pi_seq(T) = kappa * y_seq(T);
for t = T-1:-1:1
    pi_seq(t) = kappa * y_seq(t) + beta * pi_seq(t+1);
end

% monetary block

i_seq = phi_pi * pi_seq;
% i_seq = [pi_seq(2:T);0] + (phi_pi - 1) * y_seq;
% rho_tr = 0.3;
% i_seq = zeros(T,1);
% i_seq(1) = (1-rho_tr) * phi_pi * pi_seq(1);
% for t = 2:T
%     i_seq(t) = rho_tr * i_seq(t-1) + (1-rho_tr) * phi_pi * pi_seq(t);
% end

% transfers and debt

d_seq     = NaN(T,1);
t_seq     = NaN(T,1);
tau_e_seq = NaN(T,1);

% H rules

if isempty(tau_d) % this will be the case for the H rules

tau_d_aux = 1;

if H == 0

d_seq(1) = 0;
t_seq(1) = 0 + D_SS/Y_SS * (0 - pi_seq(1)) - beta * d_seq(1);
tau_e_seq(1) = t_seq(1) - tau_y * y_seq(1) - tau_d_aux * (0 + D_SS/Y_SS * (0 - pi_seq(1))) - tau_x_seq(1);

for t = 2:T
    d_seq(t) = 0;
    t_seq(t) = d_seq(t-1) + D_SS/Y_SS * (i_seq(t-1) - pi_seq(t)) - beta * d_seq(t);
    tau_e_seq(t) = t_seq(t) - tau_y * y_seq(t) - tau_d_aux * (d_seq(t-1) + D_SS/Y_SS * (i_seq(t-1) - pi_seq(t))) - tau_x_seq(t);
end

else

t_seq(1) = tau_y * y_seq(1) + tau_x_seq(1);
d_seq(1) = 1/beta * (0 - t_seq(1) + D_SS/Y_SS * (0 - pi_seq(1)));
tau_e_seq(1) = 0;

if H > 1

for t = 2:H
    t_seq(t)     = tau_y * y_seq(t) + tau_x_seq(t);
    d_seq(t)     = 1/beta * (d_seq(t-1) - t_seq(t) + D_SS/Y_SS * (i_seq(t-1)-pi_seq(t)));
    tau_e_seq(t) = 0;
end

end

for t = H+1:T
    d_seq(t) = 0;
    t_seq(t) = d_seq(t-1) + D_SS/Y_SS * (i_seq(t-1) - pi_seq(t)) - beta * d_seq(t);
    tau_e_seq(t) = t_seq(t) - tau_y * y_seq(t) - tau_d_aux * (d_seq(t-1) + D_SS/Y_SS * (i_seq(t-1) - pi_seq(t))) - tau_x_seq(t);
end

end

end

% tau_d rules

if isempty(tau_d) == 0
    
if H == 0

d_seq(1) = 0;
t_seq(1) = 0 + D_SS/Y_SS * (0 - pi_seq(1)) - beta * d_seq(1);
tau_e_seq(1) = t_seq(1) - tau_y * y_seq(1) - tau_d * (0 + D_SS/Y_SS * (0 - pi_seq(1))) - (1-tau_d) * tau_x_seq(1);

for t = 2:T
    d_seq(t) = 0;
    t_seq(t) = d_seq(t-1) + D_SS/Y_SS * (i_seq(t-1) - pi_seq(t)) - beta * d_seq(t);
    tau_e_seq(t) = t_seq(t) - tau_y * y_seq(t) - tau_d * (d_seq(t-1) + D_SS/Y_SS * (i_seq(t-1) - pi_seq(t))) - (1-tau_d) * tau_x_seq(t);
end

else

t_seq(1) = tau_y * y_seq(1) + tau_d * (0 + D_SS/Y_SS * (0 - pi_seq(1))) + (1-tau_d) * tau_x_seq(1);
d_seq(1) = 1/beta * (0 - t_seq(1) + D_SS/Y_SS * (0 - pi_seq(1)));
tau_e_seq(1) = 0;

if H > 1

for t = 2:H
    t_seq(t)     = tau_y * y_seq(t) + (1-tau_d) * tau_x_seq(t);
    d_seq(t)     = 1/beta * (d_seq(t-1) - t_seq(t) + D_SS/Y_SS * (i_seq(t-1)-pi_seq(t)));
    tau_e_seq(t) = 0;
end

end

for t = H+1:T
    t_seq(t)     = tau_y * y_seq(t) + tau_d * (d_seq(t-1) + D_SS/Y_SS * (i_seq(t-1) - pi_seq(t))) + (1-tau_d) * tau_x_seq(t);
    d_seq(t)     = 1/beta * (d_seq(t-1) - t_seq(t) + D_SS/Y_SS * (i_seq(t-1) - pi_seq(t)));
    tau_e_seq(t) = tau_d * (d_seq(t-1) + D_SS/Y_SS * (i_seq(t-1) - pi_seq(t)));
end

end

end

% consumption

c_seq    = C_y * (y_seq - t_seq) + C_i * i_seq + C_pi * pi_seq;
d_seq_HH = D_y * (y_seq - t_seq) + D_i * i_seq + D_pi * pi_seq;