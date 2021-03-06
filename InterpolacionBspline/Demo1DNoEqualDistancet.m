close all
N = 20; % upsampling rate
Noriginal = 50;
proportion = 40;



eje_x = -10:0.2:10; linspace(-10,10,Noriginal)
t = linspace(0,1,length(eje_x));
image = sinc(eje_x);

figure
subplot(3,5,1)
stem(t,image)
axis([0 1 -0.3 1])
axis tight
title('Original image without interpolation')

[A] = length(image);
[X] = linspace(0,A,N*A);
[t_prima] = linspace(0,1,N*A);
graph = [6,11];
degrees =[-3,3];
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
    subplot(3,5,graph(i))
    plot(t_prima,V_new,'LineWidth', 2)
    axis([0 1 -0.3 1])
    hold on
    stem(t,image, 'r')
    axis tight
    if SplineDeg == 3
        title_string = ['Interp. BS-deg: ' num2str(SplineDeg)];
    else
        title_string = ['Interp. oMom-deg: ' num2str(-SplineDeg)];
    end
    disp(title_string);
    title(title_string)
end

% Now we have to interpolate also the time because is not equally subsampled

% Create secuence

selection = randperm(length(eje_x));
selection = selection(1:int8((proportion * length(eje_x))/100));
selection = sort(selection);

image_sel = image(selection);
t_selection = t(selection);


subplot(3,5,2)
stem(linspace(0,1,length(image_sel)), image_sel)
axis tight
axis([0 1 -0.3 1])
title('Original Data Selection')


[A] = length(image_sel);
[X] = linspace(0,A,N*A);
[t_prima] = linspace(0,1,N*A);
graph = [7,12];
degrees =[-3,3];
timescoeffs = zeros(size(degrees)+3);
timesinterp = zeros(size(degrees)+3);
timestotal = zeros(size(degrees)+3);

for i =1:length(degrees)
    SplineDeg = degrees(i); 
    tStart=tic;
    coeffs = ConvertToInterpolationCoefficients(image_sel,SplineDeg);
    timescoeffs(i)=toc(tStart);
    V_new = X*0;
    tStart=tic;
    V_new(:) =  InterpolValue(coeffs,X',SplineDeg);
    timesinterp(i)=toc(tStart); 
    subplot(3,5,graph(i))
    plot(t_prima,V_new,'LineWidth', 2)
    hold on
    stem(linspace(0,1,length(image_sel)), image_sel, 'r')
    axis([0 1 -0.3 1])
    if SplineDeg == 3
        title_string = ['Interp. BS-deg: ' num2str(SplineDeg)];
    else
        title_string = ['Interp. oMom-deg: ' num2str(-SplineDeg)];
    end
    disp(title_string);
    title(title_string)
end

% Interpolate time

subplot(3,5,3)
stem(linspace(0,1,length(t_selection)), t_selection)
axis tight
axis([0 1 -0.3 1])
title('Original Time Selection')

[A] = length(t_selection);
[X] = linspace(0,A,N*A);
[t_prima] = linspace(0,1,N*A);
graph = [8,13];
degrees =[1,2];
timescoeffs = zeros(size(degrees)+3);
timesinterp = zeros(size(degrees)+3);
timestotal = zeros(size(degrees)+3);

for i =1:length(degrees)
    SplineDeg = degrees(i); 
    tStart=tic;
    coeffs = ConvertToInterpolationCoefficients(t_selection,SplineDeg);
    timescoeffs(i)=toc(tStart);
    V_new = X*0;
    tStart=tic;
    V_new(:) =  InterpolValue(coeffs,X',SplineDeg);
    timesinterp(i)=toc(tStart); 
    subplot(3,5,graph(i))
    plot(t_prima,V_new,'LineWidth', 2)
    axis([0 1 -0.3 1])
    hold on
    stem(linspace(0,1,length(t_selection)), t_selection, 'r')
    title_string = ['Interp. BS-deg: ' num2str(SplineDeg)];
    disp(title_string);
    title(title_string)
end

% Interpolate both

subplot(3,5,4)
stem(t_selection,image_sel)
axis tight
axis([0 1 -0.3 1])
title('Original Time Vs Data Selection')

[A] = length(t_selection);
[X] = linspace(0,A,N*A);
[t_prima] = linspace(0,1,N*A);
graph = [9,14];
degrees =[-3,3];
timescoeffs = zeros(size(degrees)+3);
timesinterp = zeros(size(degrees)+3);
timestotal = zeros(size(degrees)+3);

timeSplineDegree = 2;

for i =1:length(degrees)
    SplineDeg = degrees(i); 
    tStart=tic;
    coeffsD = ConvertToInterpolationCoefficients(image_sel,SplineDeg);
    coeffsT = ConvertToInterpolationCoefficients(t_selection,timeSplineDegree);
    timescoeffs(i)=toc(tStart);
    V_newD = X*0;
    V_newT = X*0;
    tStart=tic;
    V_newD(:) =  InterpolValue(coeffsD,X',SplineDeg);
    V_newT(:) =  InterpolValue(coeffsT,X',timeSplineDegree);
    timesinterp(i)=toc(tStart); 
    subplot(3,5,graph(i))
    plot(V_newT,V_newD,'LineWidth', 2)
    axis([0 1 -0.3 1])
    hold on
    stem(t_selection,image_sel , 'r')
    if SplineDeg == 3
        title_string = ['Interp. BS-deg: Time 2' num2str(SplineDeg)];
    else
        title_string = ['Interp. oMom-deg: Time 2' num2str(-SplineDeg)];
    end
    disp(title_string);
    title(title_string)
end

graph = [10,15];
timeSplineDegree = 1;

for i =1:length(degrees)
    SplineDeg = degrees(i); 
    tStart=tic;
    coeffsD = ConvertToInterpolationCoefficients(image_sel,SplineDeg);
    coeffsT = ConvertToInterpolationCoefficients(t_selection,timeSplineDegree);
    timescoeffs(i)=toc(tStart);
    V_newD = X*0;
    V_newT = X*0;
    tStart=tic;
    V_newD(:) =  InterpolValue(coeffsD,X',SplineDeg);
    V_newT(:) =  InterpolValue(coeffsT,X',timeSplineDegree);
    timesinterp(i)=toc(tStart); 
    subplot(3,5,graph(i))
    plot(V_newT,V_newD,'LineWidth', 2)
    axis([0 1 -0.3 1])
    hold on
    stem(t_selection,image_sel , 'r')
    if SplineDeg == 3
        title_string = ['Interp. BS-deg: Time 1' num2str(SplineDeg)];
    else
        title_string = ['Interp. oMom-deg: Time 1' num2str(-SplineDeg)];
    end
    disp(title_string);
    title(title_string)
end



timestotal = timesinterp+timescoeffs;