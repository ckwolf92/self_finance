% inflation

pi_seq = zeros(T,1);
pi_seq(T) = kappa * y_seq(T);
for t = T-1:-1:1
    pi_seq(t) = kappa * y_seq(t) + 1/(1+r_SS) * pi_seq(t+1);
end

% monetary block

i_seq = [pi_seq(2:T);0] + phi * y_seq;

% transfers and debt

d_seq     = NaN(T,1);
tau_seq   = NaN(T,1);

tau_seq(1) = tau_x_seq(1) - tau_d * (1+r_SS) * D_SS/Trans_SS * (0 + (0 - pi_seq(1)));
d_seq(1)   = (1+r_SS) * (0 + 0 - pi_seq(1)) + Trans_SS/D_SS * tau_seq(1) - tau_y * Y_SS/D_SS * y_seq(1);

for t = H+1:T
    tau_seq(t)   = tau_x_seq(t) - tau_d * (1+r_SS) * D_SS/Trans_SS * (d_seq(t-1) + (i_seq(t-1) - pi_seq(t)));
    d_seq(t)     = (1+r_SS) * (d_seq(t-1) + i_seq(t-1) - pi_seq(t)) + Trans_SS/D_SS * tau_seq(t) - tau_y * Y_SS/D_SS * y_seq(t);
end

% consumption

c_seq    = C_y * y_seq + C_tau * tau_seq + C_i * i_seq + C_pi * pi_seq;
d_seq_HH = D_y * y_seq + D_tau * tau_seq + D_i * i_seq + D_pi * pi_seq;