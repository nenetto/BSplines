% Modification to use in Matlab by:
%
%
% This code allow to interpolate with Bspline(0...9) and oMom3.
%
% E. Marinetto
% Laboratorio de Imagen Médica
% http://image.hggm.es/
% 
% Unidad de Medicina y Cirugía Experimental.
% Fundación para la Investigación Biomédica del 
% Hospital Gregorio Marañón
% C/ Doctor Esquerdo, 46
% 28007  Madrid
% Teléfono: +34 91 4265017
%  Fax: +34 91 4265108

% *----------------------------------------------------------------------------
%  *	This C program is based on the following paper:
%  *		P. Thevenaz, T. Blu, M. Unser, "Interpolation Revisited,"
%  *		IEEE Transactions on Medical Imaging,
%  *		vol. 19, no. 7, pp. 739-758, July 2000.
%  *----------------------------------------------------------------------------
%  *	Philippe Thevenaz
%  *	EPFL/STI/IMT/LIB/BM.4.137
%  *	Station 17
%  *	CH-1015 Lausanne VD
%  *----------------------------------------------------------------------------
%  *	phone (CET):	+41(21)693.51.61
%  *	fax:			+41(21)693.37.01
%  *	RFC-822:		philippe.thevenaz@epfl.ch
%  *	X-400:			/C=ch/A=400net/P=switch/O=epfl/S=thevenaz/G=philippe/
%  *	URL:			http://bigwww.epfl.ch/
%  *-----------------------------------------------------------------------
%  -----


% Makefile B-spline Thevenez

% 1) Setup the mex compiler
%
% Uncomment this line if you don't have any setup compiler first.
% You can find information and howto in: mex -help or in the web

% mex -setup

% 2) Compile necessary files
close all
clc
delete *mexw64
mex -c bspline_coeff.c
mex -c bspline_interp.c
mex ConvertToInterpolationCoefficients_c.c bspline_coeff.obj
mex InterpolValue_c.c bspline_interp.obj


% 3) Delete the unnecessary files

delete *.obj
delete *.asv

%% Test the programs 1D

% 4) Test the programs 1D

N = 4; % upsampling rate
eje_x = -10:0.4:10;
image = sinc(eje_x);

figure
subplot(4,5,3)
stem(image)
axis tight
title('Original image without interpolation')

[A] = length(image);
[X] = linspace(1,A,N*A);
graph = 6:20;
degrees =[-3,0:9];
timescoeffs = zeros(size(degrees)+3);
timesinterp = zeros(size(degrees)+3);
timestotal = zeros(size(degrees)+3);

for i =1:length(degrees)
    SplineDeg = degrees(i); 
    tStart=tic;
    coeffs = ConvertToInterpolationCoefficients(image,SplineDeg);
    timescoeffs(i)=toc(tStart);
    V_new = X*0;
    tStart=tic;
    V_new(:) =  InterpolValue(coeffs,X',SplineDeg);
    timesinterp(i)=toc(tStart); 
    subplot(4,5,graph(i))
    plot(X(:),V_new)
    axis tight
    title_string = ['Interp. BS-deg:' num2str(SplineDeg)];
    disp(title_string);
    title(title_string)

end

timestotal = timesinterp+timescoeffs;


subplot(4,5,17)
tStart=tic;
image_interp_near = interp1(double(image),X,'nearest');
timestotal(12)=toc(tStart);
plot(X(:),image_interp_near)
axis tight
title('Interpolada Vecino Más Cercano')
subplot(4,5,18)
tStart=tic;
image_interp_lin = interp1(double(image),X,'linear');
timestotal(13)=toc(tStart);
plot(X(:),image_interp_lin);
axis tight
title('Interpolada Lineal')
subplot(4,5,19)
tStart=tic;
image_interp_cubic = interp1(double(image),X,'cubic');
timestotal(14)=toc(tStart);
plot(X(:),image_interp_cubic)
axis tight
title('Interpolada Cubica')
subplot(4,5,20)
plot(degrees(2:11),timestotal(2:11),'r','LineWidth',2);hold on
scatter(3,timestotal(1),'b','LineWidth',2);
scatter([0,1,3],timestotal(12:14),'m','LineWidth',2);hold off
axis tight
title('Tiempos procesado para cada B-deg')




%% Test the programs 2D

% 4) Test the programs 2D



load clown % image to interpolate
N = 4; % upsampling rate
image = X(:,:);
clear X caption map;

figure
subplot(4,5,3)
imshow(image,[min(image(:)) max(image(:))])
title('Original image without interpolation')

[A,B] = size(image);
[X,Y] = meshgrid(linspace(1,B,N*B),linspace(1,A,N*A));
graph = 6:20;
degrees =[-3,0:9];
timescoeffs = zeros(size(degrees)+3);
timesinterp = zeros(size(degrees)+3);
timestotal = zeros(size(degrees)+3);

for i =1:length(degrees)
    SplineDeg = degrees(i); 
    tStart=tic;
    coeffs = ConvertToInterpolationCoefficients(image,SplineDeg);
    timescoeffs(i)=toc(tStart);
    V_new = X*0;
    tStart=tic;
    V_new(:) =  InterpolValue(coeffs,[Y(:),X(:)],SplineDeg);
    timesinterp(i)=toc(tStart);
    subplot(4,5,graph(i))
    imshow(V_new,[min(image(:)) max(image(:))])
    title_string = ['Interp. BS-deg:' num2str(SplineDeg)];
    disp(title_string);
    title(title_string)

end

timestotal = timesinterp+timescoeffs;


subplot(4,5,17)
tStart=tic;
image_interp_near = interp2(double(image),X,Y,'nearest');
timestotal(12)=toc(tStart);
imshow(image_interp_near,[min(image(:)) max(image(:))])
title('Interpolada Vecino Más Cercano')
subplot(4,5,18)
tStart=tic;
image_interp_lin = interp2(double(image),X,Y,'linear');
timestotal(13)=toc(tStart);
imshow(image_interp_lin,[min(image(:)) max(image(:))]);
title('Interpolada Lineal')
subplot(4,5,19)
tStart=tic;
image_interp_cubic = interp2(double(image),X,Y,'cubic');
timestotal(14)=toc(tStart);
imshow(image_interp_cubic,[min(image(:)) max(image(:))])
title('Interpolada Cubica')
subplot(4,5,20)
plot(degrees(2:11),timestotal(2:11),'r','LineWidth',2);hold on
scatter(3,timestotal(1),'b','LineWidth',2);
scatter([0,1,3],timestotal(12:14),'m','LineWidth',2);hold off
title('Tiempos procesado para cada B-deg')



%% Test the programs 3D

% 5) Test the programs 3D

CORTE = 15;

%
load mri % image to interpolate
N = 4; % upsampling rate
image = squeeze(D);
clear D siz map;

figure
subplot(4,5,3)
imshow(image(:,:,CORTE),[min(image(:)) max(image(:))])
title('Original image without interpolation')

[A,B,C] = size(image);
[X,Y,Z] = meshgrid(linspace(1,B,N*B),linspace(1,A,N*A),linspace(1,C,N*C));
graph = 6:20;
degrees =[-3,0:9];
timescoeffs = zeros(size(degrees)+3);
timesinterp = zeros(size(degrees)+3);
timestotal = zeros(size(degrees)+3);

for i =1:length(degrees)
    SplineDeg = degrees(i); 
    tStart=tic;
    coeffs = ConvertToInterpolationCoefficients(image,SplineDeg);
    timescoeffs(i)=toc(tStart);
    V_new = X*0;
    tStart=tic;
    V_new(:) =  InterpolValue(coeffs,[Y(:),X(:),Z(:)],SplineDeg);
    timesinterp(i)=toc(tStart);
    subplot(4,5,graph(i))
    imshow(V_new(:,:,CORTE),[min(image(:)) max(image(:))])
    title_string = ['Interp. BS-deg:' num2str(SplineDeg)];
    disp(title_string);
    title(title_string)

end

timestotal = timesinterp+timescoeffs;


subplot(4,5,17)
tStart=tic;
image_interp_near = interp3(double(image),X,Y,Z,'nearest');
timestotal(12)=toc(tStart);
imshow(image_interp_near(:,:,CORTE),[min(image(:)) max(image(:))])
title('Interpolada Vecino Más Cercano')
subplot(4,5,18)
tStart=tic;
image_interp_lin = interp3(double(image),X,Y,Z,'linear');
timestotal(13)=toc(tStart);
imshow(image_interp_lin(:,:,CORTE),[min(image(:)) max(image(:))]);
title('Interpolada Lineal')
subplot(4,5,19)
tStart=tic;
image_interp_cubic = interp3(double(image),X,Y,Z,'cubic');
timestotal(14)=toc(tStart);
imshow(image_interp_cubic(:,:,CORTE),[min(image(:)) max(image(:))])
title('Interpolada Cubica')
subplot(4,5,20)
plot(degrees(2:11),timestotal(2:11),'r','LineWidth',2);hold on
scatter(3,timestotal(1),'b','LineWidth',2);
scatter([0,1,3],timestotal(12:14),'m','LineWidth',2);hold off
title('Tiempos procesado para cada B-deg')