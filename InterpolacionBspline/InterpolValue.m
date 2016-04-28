function valueInterp = InterpolValue(coefficients, valuePosition, SplineDegree)

% Modificamos las variables para poder procesarlas
coefficients = double(squeeze(coefficients));
SplineDegree = round(SplineDegree);
 
if(isvector(coefficients))
    coefficients = coefficients';
end


% Check Grado Spline/oMom
if(SplineDegree>9)
    error('El grado de Bspline no está implementado');
elseif(and(SplineDegree<0,SplineDegree~=-3))
        error('El grado de oMom no está implementado');
end

[L,~] = size(valuePosition);
valueInterp = zeros(L,1);


for i=1:L
valueInterp(i) = InterpolValue_c(coefficients,valuePosition(i,:)-0.5,SplineDegree);
end

end