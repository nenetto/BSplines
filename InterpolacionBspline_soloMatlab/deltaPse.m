function y = deltaPse(n,nF,x)
y = x*0;
for k=0:n 
    y = y + ((-1)^k) * nchoosek(n,k) * pse(nF,x-k);
end
end