delete bspline_interp.obj InterpolValue_c.mexw64
mex -c bspline_interp.c
mex InterpolValue_c.c bspline_interp.obj


% Aumento de resolución usando Bsplines
clear all
close all
clc
load mri
image = squeeze(D);
image = double(image(:,:,15));


% 
%  load clown
%  image = X;

GradoSpline = 3;


%  profile on, profile clear
coeffs = ConvertToInterpolationCoefficients_c(image,GradoSpline);
%coeffs2 = SamplesToCoefficients(image,GradoSpline);
%  profile report, profile off

%
[A,B] = size(coeffs);


N = 4;
[X,Y] = meshgrid(linspace(1,B,N*B),linspace(1,A,N*A));%;,linspace(1,C,N*C));
%[X,Y,Z] = meshgrid(linspace(1,B,N*B),linspace(1,A,N*A),linspace(1,C,N*C));



V_new = X*0;


x = X(:);
y = Y(:);
% z = Z(:);
% t = T(:);


for i=1:length(X(:))
    coordenada = [y(i),x(i)];
   V_new(i) =  InterpolValue_c(coeffs,coordenada,GradoSpline);
   
end

%





%
% 2D
figure
subplot(2,3,1)
imshow(image,[min(image(:)) max(image(:))])
title('Original')
subplot(2,3,2)
imshow(V_new,[min(image(:)) max(image(:))])
title('Interpolada BSpline')
subplot(2,3,3)
image_interp_near = interp2(double(image),X,Y,'nearest');
imshow(image_interp_near,[])
title('Interpolada Vecino Más Cercano')
subplot(2,3,4)
image_interp_lin = interp2(double(image),X,Y,'linear');
imshow(image_interp_lin,[])
title('Interpolada Lineal')
subplot(2,3,5)
image_interp_cubic = interp2(double(image),X,Y,'cubic');
imshow(image_interp_cubic,[])
title('Interpolada Cubica')

%% 3D
figure
subplot(2,3,1)
imshow(image(:,:,1),[min(image(:)) max(image(:))])
title('Original')
subplot(2,3,2)
imshow(V_new(:,:,1),[min(image(:)) max(image(:))])
title('Interpolada BSpline')
subplot(2,3,3)
image_interp_near = interp2(double(image(:,:,1)),X,Y,'nearest');
imshow(image_interp_near,[])
title('Interpolada Vecino Más Cercano')
subplot(2,3,4)
image_interp_lin = interp2(double(image(:,:,1)),X,Y,'linear');
imshow(image_interp_lin,[])
title('Interpolada Lineal')
subplot(2,3,5)
image_interp_cubic = interp2(double(image(:,:,1)),X,Y,'cubic');
imshow(image_interp_cubic,[])
title('Interpolada Cubica')







