function y = lin(x)
y = zeros(size(x));
x_abs = abs(x);

y(x_abs>1) = 0;
y(x_abs<=1) = 1-x_abs(x_abs<=1);

end