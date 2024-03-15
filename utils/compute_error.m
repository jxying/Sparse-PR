function err = compute_error(x,z)
    err1  = norm(x-z)/norm(z);
    err2  = norm(x+z)/norm(z); 
    err   = min(err1,err2);    
end