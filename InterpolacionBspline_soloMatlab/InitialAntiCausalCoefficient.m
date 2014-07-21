function c = InitialAntiCausalCoefficient(samples,z)

% This initialization corresponds to mirror boundaries
DataLength = length(samples);

c = (z/(z*z -1)) * (z * samples(DataLength-1)+ samples(DataLength));

end