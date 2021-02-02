function K = objFunc2(x)
%%
% This is the objective function for filter3

m12 = x(1);
m23 = x(2);
m34 = x(3);
m45 = x(4);
    
m11 = x(5);
m22 = x(6);
m33 = x(7);
m44 = x(8);
m55 = x(9);
    
m13 = x(10);
m35 = x(11);
    
r11 = x(12);
rnn = r11;

M = [m11, m12, m13, 0, 0;
    m12, m22, m23, 0, 0;
    m13, m23, m33, m34, m35;
    0, 0, m34, m44, m45;
    0, 0, m35, m45, m55];
order = length(M);
R = zeros(order);
R(1,1) = r11;
R(order,order) = rnn;

ps = [-2.4526, -1.4324];
tz =  [-0.97388,-0.746889,-0.2485,0.4225,0.92677];
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