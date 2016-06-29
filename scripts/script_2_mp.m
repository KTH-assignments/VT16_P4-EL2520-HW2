clear
close all

sys = minphase;

wc = 0.1;
pm = pi/3;

%% 1.1
[G_num, G_denom] = tfdata(sys);

G_11 = tf(G_num{1,1}, G_denom{1,1});
G_12 = tf(G_num{1,2}, G_denom{1,2});
G_21 = tf(G_num{2,1}, G_denom{2,1});
G_22 = tf(G_num{2,2}, G_denom{2,2});

% The system transfer matrix
G = [G_11, G_12; G_21, G_22];

% The magnitude and phase of G_ii at the crossover frequency
[~, ph_1] = bode(G_11, wc);
[~, ph_2] = bode(G_22, wc);

% Compute T
T_1 = 1/wc * tan(pm - pi/2 - ph_1 * pi / 180)
T_2 = 1/wc * tan(pm - pi/2 - ph_2 * pi / 180) 

% Compute K from the Bode diagram of G_ii * F / K
s = tf('s');
l_1 = G_11 * (1 + 1 / (s * T_1));
l_2 = G_22 * (1 + 1 / (s * T_2));


[K_1_inv, ~] = bode(l_1, wc);
[K_2_inv, ~] = bode(l_2, wc);

K_1 = 1 / K_1_inv
K_2 = 1 / K_2_inv

% The final controllers
f_1 = K_1 * (1 + 1 / (s * T_1));
f_2 = K_2 * (1 + 1 / (s * T_2));

F = [f_1, 0; 0, f_2];

L = minreal(G * F);

figure
margin(L(1,1))
hold on
margin(L(2,2))


%% The sensitivity function
S = minreal(inv(eye(2) + L));

%% The complementary sensitivity function
T = minreal(inv(eye(2) + L) * L);

Gc = minreal(inv(eye(2) + G) * G);

figure
sigma(S)
hold on
sigma(T)
grid
axis([0.001 10 -60 10])

%% Step responses
sim('closedloop')
figure
plot(yout); grid