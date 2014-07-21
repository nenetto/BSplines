function coefs = ConvertToInterpolationCoefficients(samples,SplineDegree)

% Modificamos las variables para poder procesarlas
samples = double(samples);
SplineDegree = round(SplineDegree);

if(SplineDegree>9)
    error('El grado de Bspline no est� implementado');
elseif(and(SplineDegree<0,SplineDegree~=-3))
        error('El grado de oMom no est� implementado');
end
coefs = ConvertToInterpolationCoefficients_c(samples,SplineDegree);
end