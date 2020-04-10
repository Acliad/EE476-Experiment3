% Script to linearize equations in section 2 EE 424
% syms Mp Lr Jr Jp Lp g theta theta_dt theta_dt2 alpha alpha_dt alpha_dt2 ...
%     Br Bp eta_g eta_m Kg kt Vm km Rm real
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
B = [0; 0; A1(1, 4); A1(2, 4)];
C = [1 0 0 0; 0 1 0 0];
D = [0; 0];