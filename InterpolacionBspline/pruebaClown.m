% Aumento de resolución usando Bsplines
clear all
close all
clc
load mri
image = squeeze(D);
image = double(image(50:100,50:100,15));
%image = round(rand(3,3,3)*10)

% 
 load clown
 image = X;

GradoSpline = 4;


%  profile on, profile clear
coeffs = ConvertToInterpolationCoefficients_c(image,GradoSpline);
%coeffs2 = SamplesToCoefficients(image,GradoSpline);
%  profile report, profile off

%
[A,B] = size(coeffs);

N = 3;
[X,Y] = meshgrid(linspace(1,B,N*B),linspace(1,A,N*A));%;,linspace(1,C,N*C));
%

V_new = X*0;


x = X(:);
y = Y(:);
% z = Z(:);

profile on, profile clear
for i=1:length(X(:))
    coordenada = [y(i),x(i)];
   V_new(i) =  InterpolValue_c(coeffs,coordenada,GradoSpline);
end
profile report, profile off



%% 3 Dims

corte = 1;

figure
subplot(2,3,1)
imshow(image(:,:,corte),[min(image(:)) max(image(:))])
title('Original')
subplot(2,3,2)
imshow(V_new(:,:,corte),[min(image(:)) max(image(:))])
title('Interpolada BSpline')
subplot(2,3,3)
image_interp_near = interp3(double(image),X,Y,Z,'nearest');
imshow(image_interp_near(:,:,corte),[])
title('Interpolada Vecino Más Cercano')
subplot(2,3,4)
image_interp_lin = interp3(double(image),X,Y,Z,'linear');
imshow(image_interp_lin(:,:,corte),[])
title('Interpolada Lineal')
subplot(2,3,5)
image_interp_cubic = interp3(double(image),X,Y,Z,'cubic');
imshow(image_interp_cubic(:,:,corte),[])
title('Interpolada Cubica')


%%

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


%%

Diff = abs(image_interp_lin-V_new);
figure
subplot(1,3,1)
imshow(V_new,[0 200])
title('Interpolada Eu')
subplot(1,3,2)
imshow(image_interp_lin,[0 200])
title('Interpolada Lineal')
subplot(1,3,3)
imshow(Diff,[0 200]);


%%
DiffA = abs(V_new-image_interp_near);
DiffB = abs(V_new-image_interp_lin);
DiffC = abs(V_new-image_interp_cubic);
close all
subplot(1,3,1)
imshow(DiffA,[])
title('Eu-near')
subplot(1,3,2)
imshow(DiffB,[])
title('Eu-lin')
subplot(1,3,3)
imshow(DiffC,[])
title('Eu-cubic')





