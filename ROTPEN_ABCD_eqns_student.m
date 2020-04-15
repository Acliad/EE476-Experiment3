syms Vm theta theta_dt theta_dt2 alpha alpha_dt alpha_dt2 real

tau = eta_g*Kg*eta_m*kt*(Vm - Kg * km * theta_dt)/Rm;

f1 = (Mp*Lr^2 + .25*Mp*Lp^2 - .25*Mp*Lp^2*cos(alpha)^2+Jr)*theta_dt2 ...
    - (.5*Mp*Lp*Lr*cos(alpha))*alpha_dt2 ...
    + (.5*Mp*Lp^2*sin(alpha)*cos(alpha)) * theta_dt * alpha_dt ...
    + (.5*Mp*Lp*Lr*sin(alpha))*alpha_dt^2;

f2 = -0.5*Mp*Lp*Lr*cos(alpha)*theta_dt2 ...
    + (Jp + 0.25 * Mp *Lp^2) * alpha_dt2 ...
    - 0.25*Mp*Lp^2*cos(alpha)*sin(alpha)*theta_dt^2 ...
    - 0.5 * Mp *Lp * g * sin(alpha);

state_vars_all = [theta theta_dt theta_dt2 alpha alpha_dt alpha_dt2];
state_vars = [theta alpha theta_dt alpha_dt ];

j = jacobian([f1, f2], state_vars_all);

j_lin = subs(j, [theta, theta_dt, alpha, alpha_dt], [0 0 0 0]);

f1_lin = j_lin(1, :) * state_vars_all' + Br * theta_dt == tau;
f2_lin = j_lin(2, :) * state_vars_all' + Bp * alpha_dt == 0;

D = [j_lin(1, 3) j_lin(1, 6); j_lin(2, 3) j_lin(2, 6)];
C = [Br 0; j_lin(2, 2) Bp];
g = [j_lin(1, 1) * theta + j_lin(2, 1) * alpha; ...
    j_lin(1, 4) * theta + j_lin(2, 4) * alpha];

S = solve([f1_lin f2_lin], [theta_dt2 alpha_dt2]);

theta_dt2_lin = S.theta_dt2 == theta_dt2;
alpha_dt2_lin = S.alpha_dt2 == alpha_dt2;

[A1, b] = equationsToMatrix([theta_dt2_lin, alpha_dt2_lin], ...
                            [state_vars, Vm]);

%  | theta  |
%  | alpha  |
%  |theta_dt|
%  |alpha_dt|
A = [0 0 1 0; 0 0 0 1; A1(1, 1:4); A1(2, 1:4)];
% TODO: Solve correctly
% B = [0; 0; A1(1, 4); A1(2, 4)];
B = [0; 0; 83.4659; 80.3162];
C = [1 0 0 0; 0 1 0 0];
D = [0; 0];

% Convert to numerical values

A = double(A);
B = double(B);
sys = ss(A, B, C, D);

%%

char_coef = poly(A); % Characteristic equation coefficients 

length_n = length(char_coef(2:end));

A_companion = zeros(length_n, length_n);

% Set up 1s
A_companion(1:end-1, 2:end) = eye(length_n-1, length_n-1);
A_companion(end, :) = -char_coef(end:-1:2);

% Set up B vector
B_companion = zeros(length_n, 1);
B_companion(end) = 1;

zeta_desired = 0.7;
wn_desired = 4; % rad/s

p1_desired = -wn_desired*exp(1j*acos(zeta_desired));
p2_desired = conj(p1_desired);
p3_desired = -30;
p4_desired = -40;

% Characteristic equation coefficients 
char_desired_coef = conv([1 -p1_desired], [1 -p2_desired]);
char_desired_coef = conv(char_desired_coef, [1 -p3_desired]);
char_desired_coef = conv(char_desired_coef, [1 -p4_desired]);

K_bar = char_desired_coef(2:end) - char_coef(2:end);

A_companion_desired = A_companion;
A_companion_desired(end, :) =  A_companion(end, :) - flip(K_bar);

T = ctrb(A, B);

T_bar = [
         B_companion, A_companion * B_companion, ...
         A_companion^2 * B_companion, A_companion^3 * B_companion
         ];

W = T/(T_bar);

K = flip(K_bar)/W;