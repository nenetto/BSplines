function image_out = InitialAntiCausalCoefficient_vector(image,z)

image_out = image;

cte = z/(z*z-1);
image_out(end,:,:,:) = cte * image(end,:,:,:) + cte*z*image(end-1,:,:,:) ;
end