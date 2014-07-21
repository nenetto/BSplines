function image_out = InitialCausalCoefficient_vector(image,z)

[Size_data(1) Size_data(2) Size_data(3) Size_data(4)] = size(image);
image_out = image;

% This initilization corresponds to mirror boundaries
TruncatedSum = 0;
if(eps > 0)
    Horizon = ceil(log(eps)/log(abs(z)));
    TruncatedSum = Horizon < Size_data(1);
end

if(TruncatedSum)
    % Accelerated loop
    if(Horizon>0)
        zn = repmat(z.^(0:Horizon)',[1 Size_data(2) Size_data(3) Size_data(4)]);
        Sum = image(1:(Horizon+1),:,:,:);
        Sum = zn.*Sum;
        Sum = sum(Sum,1);
        image(1,:,:,:) = Sum;
    end
    image_out = image;
else
    % Full loop NOT FIXED
    
    zn = repmat(z.^(0:(Size_data(1)-1))',[1 Size_data(2) Size_data(3) Size_data(4)]);
    image(1,:,:,:) = image(1,:,:,:)+z^(Size_data(1)-1)*image(end,:,:,:);
    Sum = image(1:end,:,:,:);
    z2n = z^(2*Size_data(1)-2)*(1./zn);
    Sum = (zn+z2n).*Sum;
    Sum = sum(Sum,1);
    image_out(1,:,:,:) = Sum/(1-z^(2*Size_data(1)-2));
end

end