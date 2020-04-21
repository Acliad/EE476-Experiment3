% IMPORTANT: Make sure you run setup_rotpen.m first. You need the (A,B)
% state-space matrices.
%
% Find poles of closed-loop system. 
% Verify that they are the same as the desired poles.

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