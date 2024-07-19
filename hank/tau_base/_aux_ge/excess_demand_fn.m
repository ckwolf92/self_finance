function excess_demand = excess_demand_fn(guess_seq);

%% GLOBAL VARIABLES

global beta gamma eta alpha tau_y r_SS D_SS Y_SS Trans_SS D_i D_pi D_y D_tau C_i C_pi C_y C_tau kappa ...
    phi tau_d nM_indic H tau_x_seq ...
    T

%% COLLECT INPUTS

y_seq  = guess_seq(1:T,1);

%% GET OUTCOMES

get_aggregates

%% CHECK ACCURACY

excess_demand = d_seq_HH - d_seq;