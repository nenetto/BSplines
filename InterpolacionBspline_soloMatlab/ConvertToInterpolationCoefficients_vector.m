function coefficients  = ConvertToInterpolationCoefficients_vector(image,poles)

coefficients = image;
coefficients = squeeze(coefficients);
N = ndims(coefficients);
index = (1:N)';
Size_img = size(coefficients);
NbPoles = numel(poles);

for i=1:N
    image_aux = permute(coefficients,circshift(index,i)); 
    if(Size_img(i) ~= 1) % Special Case required by mirror boundaries
 
        % Overall Gain
        lambda = cumprod((1-poles).*(1-1./poles));
        for z=1:NbPoles %FILTERING
            % Initial Cassual coefficient
            image_aux = InitialCausalCoefficient_vector(image_aux,poles(z));
            % Causal Recursion           
            image_aux = filter([1 poles(z)],1,image_aux);   
            % Anticausal initialization
            image_aux = InitialAntiCausalCoefficient_vector(image_aux,poles(z));
            % Anticausal Recursion
            image_aux = flipdim(filter([1 poles(z)],1,flipdim(image_aux,1)),1); 
        end
        % Apply Gain
        image_aux = lambda * image_aux;
    end
    coefficients = permute(image_aux,circshift(index,i));  
end

coefficients = reshape(coefficients,size(image));

end