function INDEX = computeInterpolationIndex(SplineDeg,valuePosition)


W = SplineDeg+1;

% OMoms3
if(SplineDeg == -1)
    W = 4;
end

Dmax = W/2-eps;
dimensions = length(valuePosition);
INDEX = zeros(dimensions,W);



for k = 1:dimensions

    i = ceil(valuePosition(k)-Dmax);
    INDEX(k,:) = i:(i+W-1);

end
    


end