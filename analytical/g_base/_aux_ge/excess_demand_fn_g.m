function excess_demand = excess_demand_fn_g(guess_seq);

%% GLOBAL VARIABLES

global tau_y beta r_SS D_SS Y_SS D_i D_pi D_y C_i C_pi C_y kappa ...
    phi tau_d nM_indic H g_seq ...
    T

%% COLLECT INPUTS

y_seq  = guess_seq(1:T,1);

%% GET OUTCOMES

get_aggregates_g

%% CHECK ACCURACY

% excess_demand = c_seq + g_seq - y_seq;
excess_demand = d_seq_HH - d_seq;