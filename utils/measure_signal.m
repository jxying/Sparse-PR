function [y_abs, A] = measure_signal(m,z)
n     = length(z);

%% signal measurements
A     = randn(m,n);   % sensing matrix
y     = A*z;        
y_abs = abs(y);   % measurements

end