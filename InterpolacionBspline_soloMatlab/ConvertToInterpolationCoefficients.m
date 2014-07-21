function coefs = ConvertToInterpolationCoefficients(samples,poles,tol)


lambda = 1.0;
DataLength = length(samples);
NbPoles = numel(poles);


coefs = zeros(size(samples));


if(DataLength == 1) % Special Case required by mirror boundaries
    coefs = samples;
else
    % Compute the overall gain
   for k=1:NbPoles
       lambda = lambda * (1-poles(k))*(1-1/poles(k));
   end

   %loop over all poles
   for k=1:NbPoles
      % Causal initialization
      coefs(1) = InitialCausalCoefficient(samples,poles(k),tol);
      % Causal Recursion
      for n=2:DataLength
          coefs(n) = samples(n) + poles(k)*coefs(n-1);
      end
      % Anticausal initialization
      coefs(DataLength) = InitialAntiCausalCoefficient(samples,poles(k));
      for n=(DataLength-1):-1:1
          coefs(n) = poles(k) * (coefs(n+1)-coefs(n));
      end
   end
   
      % Apply the gain
   for n=1:DataLength
       coefs(n) = coefs(n) * lambda;
   end

end
end