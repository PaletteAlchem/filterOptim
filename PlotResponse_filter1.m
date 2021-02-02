function [w, S21,S11] = PlotResponse_filter1(x)
m12 = x(1);
m23 = x(2);
m34 = x(3);

m14 = x(4);
r11 = x(5);
rnn = r11;

M = [0, m12, 0, m14;
     m12, 0, m23, 0;
     0, m23, 0, m34;
     m14, 0, m34, 0];
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