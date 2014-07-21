% Demo en 3D con Bspline 3 y oMoms3

CORTE = 15;
load mri % image to interpolate
N = 4; % upsampling rate
image = squeeze(D);
clear D siz map;

figure
imshow(image(:,:,CORTE),[min(image(:)) max(image(:))])
title('Original image without interpolation')


CORTE2 = CORTE*4;


SplineDeg = 3; 
tStart=tic;
coeffs = ConvertToInterpolationCoefficients(image,SplineDeg);
timescoeffs(1)=toc(tStart);
V_new = X*0;
tStart=tic;
V_new(:) =  InterpolValue(coeffs,[Y(:),X(:),Z(:)],SplineDeg);
timesinterp(1)=toc(tStart);
figure
imshow(V_new(:,:,CORTE2),[min(image(:)) max(image(:))])
title_string = ['Interp. BSpl 3'];
disp(title_string);
title(title_string)

SplineDeg = -3; 
tStart=tic;
coeffs = ConvertToInterpolationCoefficients(image,SplineDeg);
timescoeffs(2)=toc(tStart);
V_new = X*0;
tStart=tic;
V_new(:) =  InterpolValue(coeffs,[Y(:),X(:),Z(:)],SplineDeg);
timesinterp(2)=toc(tStart);
figure
imshow(V_new(:,:,CORTE2),[min(image(:)) max(image(:))])
title_string = ['Interp. oMom 3'];
disp(title_string);
title(title_string)


disp(['Tiempo Coeficientes BSpl 3:',num2str(timescoeffs(1))]);
disp(['Tiempo Coeficientes oMom 3:',num2str(timescoeffs(2))]);
disp(['Tiempo Interpolatio BSpl 3:',num2str(timesinterp(1))]);
disp(['Tiempo Interpolatio oMom 3:',num2str(timesinterp(2))]);


