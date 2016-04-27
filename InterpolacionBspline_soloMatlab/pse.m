function y = pse(n,x)
    n = double(n);

    if(n>=0)
        y = ((x.^n).*sign(x))./(2*factorial(n));
    else
        y = 0*x;
    
    end
    
end