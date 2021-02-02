function [w, S21,S11] = PlotResponse_filter3(x)
%% Calculate the responses (S-parameters) from optimized coupling matrix;
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


w = linspace(-5,5,500);
for index = 1:length(w)
    A = A_matrix(w(index), M, R);
    [S21(index), S11(index)] = M2SP(A,R);
end
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