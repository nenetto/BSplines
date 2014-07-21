function Pesos = computeInterpolationWeights(SplineDeg,INDEX,valuePosition)

[dimensions,W] = size(INDEX); 
Pesos = zeros(dimensions,W);

for k=1:dimensions
%     for i=1:W
%         w = valuePosition(k) - INDEX(k,i);
%         %Pesos(k,i) = Bspline(SplineDeg,w);
%         
%     end
    Pesos(k,:) = BsplineWeigth(SplineDeg,INDEX(k,:),valuePosition(k));
end
             
        
end


