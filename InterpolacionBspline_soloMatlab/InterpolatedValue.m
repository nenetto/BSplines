function valueInterp = InterpolatedValue(coefficients, valuePosition, SplineDegree)

[Nlines,Ncolumns,Nslices,Nframes] = size(coefficients);
size_img = [Nlines,Ncolumns,Nslices,Nframes];



valueInterp = 0;
index = computeInterpolationIndex(SplineDegree,valuePosition);
pesos = computeInterpolationWeights(SplineDegree,index,valuePosition);

[dimensions,W] = size(index); 

for a = 1:dimensions  
    L = size_img(a); 
    index(a,1:W) = correctIndex(index(a,1:W),L);
end

switch(dimensions)
    case 1
        for i=1:W
            valueInterp = valueInterp + sum(coefficients(index(i))*pesos(i));
        end
    case 2
        for i=1:W
            for j=1:W    
                valueInterp = valueInterp + sum(coefficients(index(1,i),index(2,j))*pesos(1,i)*pesos(2,j));
            end
        end
    case 3
        for i=1:W
            for j=1:W
                for k=1:W
                    valueInterp = valueInterp + sum(coefficients(index(1,i),index(2,j),index(3,k))*pesos(1,i)*pesos(2,j)*pesos(3,k));
                end
            end
        end
    case 4
        for i=1:W
            for j=1:W
                for k=1:W
                    for t=1:W
                    valueInterp = valueInterp + sum(coefficients(index(1,i),index(2,j),index(3,k),index(4,t))*pesos(1,i)*pesos(2,j)*pesos(3,k)*pesos(4,t));
                    end
                end
            end
        end

end




end