% Andrei Leonard Nicusan, 2018

function [ I1 ] = IntegralNC(f,a,b,m,n)
% This is a Matlab function that numerically calculates the definite
% integral of a function using the Newton-Cotes Quadrature Formula, for
% any given number of divisions, and any given degree.
    % f = function
    % a,b = interval
    % m = number of divisions
    % n = order (3 = Newton)
I1 = 0;
for l = 0:(m - 1)
    y0 = a + l * (b - a) / m;
    y1 = a + (l + 1) * (b - a) / m;
    interval = (y1 - y0) / n;
    Idiv = 0;
    for i = 0:n
        syms t
        Fsyms = 1;
        for j=0:n
            Fsyms = Fsyms .* (t - j);
        end
        Fsyms = Fsyms ./ (t - i);
        Ffun = matlabFunction(Fsyms);
        w = 1/n*(-1)^(n-i)/factorial(n-i)/factorial(i)*integral(Ffun,0,n);
        Idiv = Idiv + feval(f, y0 + i * interval) * w;
    end
    Idiv = Idiv * (y1 - y0);
    I1 = I1 + Idiv;
end
end