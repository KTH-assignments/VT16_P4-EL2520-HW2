clear
close all

sys = minphase
% sys = nonminphase

%% 1.1
[G_num, G_denom] = tfdata(sys);

G_11 = tf(G_num{1,1}, G_denom{1,1});
G_12 = tf(G_num{1,2}, G_denom{1,2});
G_21 = tf(G_num{2,1}, G_denom{2,1});
G_22 = tf(G_num{2,2}, G_denom{2,2});

% The system transfer matrix
G = [G_11, G_12; G_21, G_22];

% Poles of the individual elements
p_11 = pole(G_11)
p_12 = pole(G_12)
p_21 = pole(G_21)
p_22 = pole(G_22)

%% 1.2

% Poles of G
e = eig(sys.a)

% Zeros of G
G_z = G_11 * G_22 - G_12 * G_21;

tzero(G_z)

%% 1.3

% Plot min and max singular values
figure
sigma(sys.a, sys.b, sys.c, sys.d)

% Singular values at frequency 0
[U,S,V] = svd(evalfr(G,0))

%% 1.4

RGA = evalfr(G,0) .* inv(evalfr(G,0))'

%% 1.5

figure
step(G)
