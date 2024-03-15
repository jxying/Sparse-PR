function z = generate_true_signal(n,s)
        supp    = randperm(n,s);   % generate support of the true signal
        z       = zeros(n,1);
        z(supp) = randn(s,1);   % generate s-sparse true signal
        z       = z/norm(z);          
end