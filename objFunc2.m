function K = objFunc3(x)
%%
% This is the objective function for filter2;

m12 = x(1);
m23 = x(2);
m34 = x(3);
m45 = x(4);
m56 = x(5);
    
m16 = x(6);
m25 = x(7);
    
r11 = x(8);
rnn = r11;

M = [0, m12, 0, 0, 0, m16;
     m12, 0, m23, 0, m25, 0;
     0, m23, 0, m34, 0, 0;
     0, 0, m34, 0, m45, 0;
     0, m25, 0, m45, 0, m56;
     m16, 0, 0, 0, m56, 0 ];
order = length(M);
R = zeros(order);
R(1,1) = r11;
R(order,order) = rnn;

ps = [-3.40000000000000, -1.80000000000000, 1.80000000000000, 3.40000000000000];
tz =  [-0.970705180757279, -0.733622338072634, -0.277207807158598,0.277207807158597, 0.733622338072634, 0.970705180757278];

RL = 20;
epsilon = (10^(RL/10)-1)^(-1/2);

for i = 1:length(tz)
    A = A_matrix(tz(i), M, R);
    [~,s11] = M2SP(A, R);
    K_part1_temp(i) = (abs(s11)^2);
end

for i = 1:length(ps)
    A = A_matrix(ps(i), M, R);
    [s21,~] = M2SP(A, R);
    K_part2_temp(i) = (abs(s21)^2);
end

K_part1 = sum(K_part1_temp);

K_part2 = sum(K_part2_temp);

A = A_matrix(+1, M, R);
[~,s11] = M2SP(A, R);
K_part3 = (abs(s11) - epsilon/sqrt(1+epsilon^2))^2;


A = A_matrix(-1, M, R);
[~,s11] = M2SP(A, R);
K_part4 = (abs(s11) - epsilon/sqrt(1+epsilon^2))^2;

K = K_part1 + K_part2 + K_part3 + K_part4;


end

function A = A_matrix(w, M, R)
ONE = eye(length(M));
A = ONE * w - 1j * R + M;
end

function [s21,s11] = M2SP(A,R)
r11 = R(1,1);
rnn = R(end,end);
A_inv = pinv(A);
s21 = A_inv(length(A), 1) * -2j * sqrt(r11*rnn);
s11 = 1 + 2j * r11 * A_inv(1,1);
end