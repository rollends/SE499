%% Construct symbolic quintic spline.
clear all;
close all;

A = sym('A', [1 6]);
B = sym('B', [1 6]);
x = sym('x', [1 2]);
t = sym('t');
tp = t .^ [5 4 3 2 1 0];

assume(A, 'real');
assume(B, 'real');
assume(t, 'real');
assumeAlso(t >= 0);
assumeAlso(t <= 1);
assume(x, 'real');

A(6) = A(6) - x(1);
B(6) = B(6) - x(2);

orderB = length(B) - 1;
orderA = length(A) - 1;
totalOrder = orderA + orderB;
baseA = [A zeros(1, totalOrder - length(A))];
baseB = [B zeros(1, totalOrder - length(B))];

Syl = sym(zeros(orderA + orderB));
for i = 1:orderB
    Syl(i,:) = circshift(baseA, i-1);
end
for i = 1:orderA
    Syl(orderB+i,:) = circshift(baseB, i-1);
end

DetSyl = det(Syl);

%% Generate Matlab Functions
matlabFunction(collect(DetSyl, x(1)), 'Optimize', false, 'File', 'levelquintic');
matlabFunction(collect(jacobian(DetSyl, x), x(1)), 'Optimize', false, 'File', 'jacobianquintic');
matlabFunction(collect(hessian(DetSyl, x), x(1)), 'Optimize', false, 'File', 'hessianquintic');

%% Generate C Functions (That we will compile)
ccode(collect(DetSyl, x(1)), 'File', 'levelquintic');
ccode(collect(jacobian(DetSyl, x), x(1)), 'File', 'jacobianquintic');
ccode(collect(hessian(DetSyl, x), x(1)), 'File', 'hessianquintic');


%% Construct symbolic quartic spline.
A = sym('A', [1 5]);
B = sym('B', [1 5]);
x = sym('x', [1 2]);
t = sym('t');
tp = t .^ [4 3 2 1 0];

A(5) = A(5) - x(1);
B(5) = B(5) - x(2);

orderB = length(B) - 1;
orderA = length(A) - 1;
totalOrder = orderA + orderB;
baseA = [A zeros(1, totalOrder - length(A))];
baseB = [B zeros(1, totalOrder - length(B))];

Syl = sym(zeros(orderA + orderB));
for i = 1:orderB
    Syl(i,:) = circshift(baseA, i-1);
end
for i = 1:orderA
    Syl(orderB+i,:) = circshift(baseB, i-1);
end

DetSyl = det(Syl);
matlabFunction(horner(DetSyl), 'Optimize', false, 'File', 'levelquartic')