function coefficients = SamplesToCoefficients(imagen,SplineDegree)

poles = [];
imagen = double(imagen);

switch(SplineDegree)
    case 0    
        coefficients = imagen;
        return;
        
    case 1       
        coefficients = imagen;
        return;
        
    case 2
        
        poles(1) = sqrt(8)-3;
        
    case 3
        
        poles(1) = sqrt(3)-2;
        
    case 4

        poles(1) = sqrt(664.0 - sqrt(438976.0)) + sqrt(304.0) - 19.0;
        poles(2) = sqrt(664.0 + sqrt(438976.0)) - sqrt(304.0) - 19.0;

    case 5

        poles(1) = sqrt(135.0 / 2.0 - sqrt(17745.0 / 4.0)) + sqrt(105.0 / 4.0)- 13.0 / 2.0;
        poles(2) = sqrt(135.0 / 2.0 + sqrt(17745.0 / 4.0)) - sqrt(105.0 / 4.0)- 13.0 / 2.0;

    case 6

        poles(1) = -0.48829458930304475513011803888378906211227916123938;
        poles(2) = -0.081679271076237512597937765737059080653379610398148;
        poles(3) = -0.0014141518083258177510872439765585925278641690553467;

    case 7

        poles(1) = -0.53528043079643816554240378168164607183392315234269;
        poles(2) = -0.12255461519232669051527226435935734360548654942730;
        poles(3) = -0.0091486948096082769285930216516478534156925639545994;

    case 8

        poles(1) = -0.57468690924876543053013930412874542429066157804125;
        poles(2) = -0.16303526929728093524055189686073705223476814550830;
        poles(3) = -0.023632294694844850023403919296361320612665920854629;
        poles(4) = -0.00015382131064169091173935253018402160762964054070043;

    case 9
        
        poles(1) = -0.60799738916862577900772082395428976943963471853991;
        poles(2) = -0.20175052019315323879606468505597043468089886575747;
        poles(3) = -0.043222608540481752133321142979429688265852380231497;
        poles(4) = -0.0021213069031808184203048965578486234220548560988624;  
        
    case -3
        
        poles(1) = (1/8)*(sqrt(105)-13);
        
    otherwise
        disp('Invalid Degree');
        
end

% convert the image samples into interpolation coefficients 
	% in-place separable process, along x 
        
        coefficients = ConvertToInterpolationCoefficients_vector(imagen,poles);
 

end