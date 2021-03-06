function valueInterp = InterpolValue(coefficients, valuePosition, SplineDegree)

% Modificamos las variables para poder procesarlas
coefficients = double(squeeze(coefficients));
valuePosition = double(squeeze(valuePosition));
SplineDegree = round(SplineDegree);
 
if(isvector(coefficients))
    coefficients = coefficients';
end

if(isvector(valuePosition))
    valuePosition = valuePosition';
end


% Check Grado Spline/oMom
if(SplineDegree>9)
    error('El grado de Bspline no est� implementado');
elseif(and(SplineDegree<0,SplineDegree~=-3))
        error('El grado de oMom no est� implementado');
end

% Check if column

[L,~] = size(valuePosition);
valueInterp = zeros(L,1);


for i=1:L
valueInterp(i) = InterpolValue_c(coefficients,valuePosition(i,:)-0.5,SplineDegree);
end

end