/*****************************************************************************
 * *
 *	Date: February 20, 2013
 *
 *	Some modifications over previous version were made to accomplish a
 *	N-Dimensional interpolation. Currently, this software produce Bspline
 *	Coefficients and oMom3 coefficients for interpolation up to 4 dimensions.
 *
 *	--------------------------------------------------------------------------
		E. Marinetto
		Laboratorio de Imagen Médica
		WEB del LIM

		Unidad de Medicina y Cirugía Experimental.
		Fundación para la Investigación Biomédica del
		Hospital Gregorio Marañón
		C/ Doctor Esquerdo, 46
		28007  Madrid
		Teléfono: +34 91 4265017
		 Fax: +34 91 4265108
 *
 *
 *
 *	Date: January 5, 2009
 *----------------------------------------------------------------------------
 *	This C program is based on the following paper:
 *		P. Thevenaz, T. Blu, M. Unser, "Interpolation Revisited,"
 *		IEEE Transactions on Medical Imaging,
 *		vol. 19, no. 7, pp. 739-758, July 2000.
 *----------------------------------------------------------------------------
 *	Philippe Thevenaz
 *	EPFL/STI/IMT/LIB/BM.4.137
 *	Station 17
 *	CH-1015 Lausanne VD
 *----------------------------------------------------------------------------
 *	phone (CET):	+41(21)693.51.61
 *	fax:			+41(21)693.37.01
 *	RFC-822:		philippe.thevenaz@epfl.ch
 *	X-400:			/C=ch/A=400net/P=switch/O=epfl/S=thevenaz/G=philippe/
 *	URL:			http://bigwww.epfl.ch/
 *----------------------------------------------------------------------------
 *----------------------------------------------------------------------------
 *----------------------------------------------------------------------------
 ****************************************************************************/


/*****************************************************************************
 *	Includes
 ****************************************************************************/
#include	"bspline_coeff.h"
 


/*****************************************************************************
 *	Definition of static procedures
 ****************************************************************************/
/*--------------------------------------------------------------------------*/
static void		ConvertToInterpolationCoefficients
				(
					double	c[],		/* input samples --> output coefficients */
					long	DataLength,	/* number of samples or coefficients */
					double	z[],		/* poles */
					long	NbPoles,	/* number of poles */
					double	Tolerance	/* admissible relative error */
				)

{ /* begin ConvertToInterpolationCoefficients */

	double	Lambda = 1.0;
	long	n, k;

	/* special case required by mirror boundaries */
	if (DataLength == 1L) {
		return;
	}
	/* compute the overall gain */
	for (k = 0L; k < NbPoles; k++) {
		Lambda = Lambda * (1.0 - z[k]) * (1.0 - 1.0 / z[k]);
        	}
	/* apply the gain */
	for (n = 0L; n < DataLength; n++) {
		c[n] *= Lambda;
	}
    
	/* loop over all poles */
	for (k = 0L; k < NbPoles; k++) {
		/* causal initialization */
		c[0] = InitialCausalCoefficient(c, DataLength, z[k], Tolerance);
        		/* causal recursion */
		for (n = 1L; n < DataLength; n++) {
			c[n] += z[k] * c[n - 1L];
            		}
		/* anticausal initialization */
		c[DataLength - 1L] = InitialAntiCausalCoefficient(c, DataLength, z[k]);
		/* anticausal recursion */
		for (n = DataLength - 2L; 0 <= n; n--) {
			c[n] = z[k] * (c[n + 1L] - c[n]);
		}
	}
} /* end ConvertToInterpolationCoefficients */

/*--------------------------------------------------------------------------*/
static double	InitialCausalCoefficient
				(
					double	c[],		/* coefficients */
					long	DataLength,	/* number of coefficients */
					double	z,			/* actual pole */
					double	Tolerance	/* admissible relative error */
				)

{ /* begin InitialCausalCoefficient */

	double	Sum, zn, z2n, iz;
	long	n, Horizon;

	/* this initialization corresponds to mirror boundaries */
	Horizon = DataLength;
	if (Tolerance > 0.0) {
		Horizon = (long)ceil(log(Tolerance) / log(fabs(z)));
	}
	if (Horizon < DataLength) {
		/* accelerated loop */
		zn = z;
		Sum = c[0];
		for (n = 1L; n < Horizon; n++) {
			Sum += zn * c[n];
			zn *= z;
		}
		return(Sum);
	}
	else {
		/* full loop */
		zn = z;
		iz = 1.0 / z;
		z2n = pow(z, (double)(DataLength - 1L));
		Sum = c[0] + z2n * c[DataLength - 1L];
		z2n *= z2n * iz;
		for (n = 1L; n <= DataLength - 2L; n++) {
			Sum += (zn + z2n) * c[n];
			zn *= z;
			z2n *= iz;
		}
		return(Sum / (1.0 - zn * zn));
	}
} /* end InitialCausalCoefficient */

/*--------------------------------------------------------------------------*/
static void		GetColumn
				(
					float	*Image,		/* input image array */
					long	Width,		/* width of the image */
					long	x,			/* x coordinate of the selected line */
					double	Line[],		/* output linear array */
					long	Height		/* length of the line */
				)

{ /* begin GetColumn */

	long	y;

	Image = Image + (ptrdiff_t)x;
	for (y = 0L; y < Height; y++) {
		Line[y] = (double)*Image;
		Image += (ptrdiff_t)Width;
	}
} /* end GetColumn */

/*--------------------------------------------------------------------------*/
static void		GetRow
				(
					float	*Image,		/* input image array */
					long	y,			/* y coordinate of the selected line */
					double	Line[],		/* output linear array */
					long	Width		/* length of the line */
				)

{ /* begin GetRow */

	long	x;

	Image = Image + (ptrdiff_t)(y * Width);
	for (x = 0L; x < Width; x++) {
		Line[x] = (double)*Image++;
	}
} /* end GetRow */

/*--------------------------------------------------------------------------*/
static double	InitialAntiCausalCoefficient
				(
					double	c[],		/* coefficients */
					long	DataLength,	/* number of samples or coefficients */
					double	z			/* actual pole */
				)

{ /* begin InitialAntiCausalCoefficient */

	/* this initialization corresponds to mirror boundaries */
	return((z / (z * z - 1.0)) * (z * c[DataLength - 2L] + c[DataLength - 1L]));
} /* end InitialAntiCausalCoefficient */

/*--------------------------------------------------------------------------*/
static void		PutColumn
				(
					float	*Image,		/* output image array */
					long	Width,		/* width of the image */
					long	x,			/* x coordinate of the selected line */
					double	Line[],		/* input linear array */
					long	Height		/* length of the line and height of the image */
				)

{ /* begin PutColumn */

	long	y;

	Image = Image + (ptrdiff_t)x;
	for (y = 0L; y < Height; y++) {
		*Image = (float)Line[y];
		Image += (ptrdiff_t)Width;
	}
} /* end PutColumn */

/*--------------------------------------------------------------------------*/
static void		PutRow
				(
					float	*Image,		/* output image array */
					long	y,			/* y coordinate of the selected line */
					double	Line[],		/* input linear array */
					long	Width		/* length of the line and width of the image */
				)

{ /* begin PutRow */

	long	x;

	Image = Image + (ptrdiff_t)(y * Width);
	for (x = 0L; x < Width; x++) {
		*Image++ = (float)Line[x];
	}
} /* end PutRow */

/*--------------------------------------------------------------------------*/

static double*	GetMetaRow
				(
					double	*Image,			/* input image array */
					long	*Size,			/* array with the dimensions of the image */
					long	*coordinate,	/* Meta coordinate of the selected line */
					long	direction,		/* Direction of the row */
					long	ImageDim		/* Number of dimensions of the image */
				)
{ /* begin GetMetaRow */


	long i,j;
	// Creation of the line
	double* Line = (double*) calloc(Size[direction],(int)sizeof(double));
	long jump_aux;
	long jump = 0L;
	long salto_dir = 1L;

	for(i=0L;i<ImageDim;i++){
		if(i==direction){
			jump_aux = 0L;
		}else{
			jump_aux = coordinate[i];
			for(j=i-1;j>=0L;j--){
				jump_aux = jump_aux * Size[j];
			}
		}
		jump = jump +jump_aux;
	}

	for(i=0L;i<direction;i++){
		salto_dir = salto_dir * Size[i];
	}


	for(i=0L;i<Size[direction];i++){
		Line[i] = Image[i*salto_dir+jump];
	}

	return Line;

} /* end GetMetaRow */
/*--------------------------------------------------------------------------*/
static void 	PutMetaRow
				(
					double	*Image,			/* input image array */
					long	*Size,			/* array with the dimensions of the image */
					long	*coordinate,	/* Meta coordinate of the selected line */
					long	direction,		/* Direction of the row */
					long	ImageDim,		/* Number of dimensions of the image */
					double  Line[]			/* Line to put in*/
				){

	long i,j;
	// Creation of the line
	long jump_aux;
	long jump = 0L;
	long salto_dir = 1L;

	for(i=0L;i<ImageDim;i++){
		if(i==direction){
			jump_aux = 0L;
		}else{
			jump_aux = coordinate[i];
			for(j=i-1;j>=0L;j--){
				jump_aux = jump_aux * Size[j];
			}
		}
		jump = jump +jump_aux;
	}

	for(i=0L;i<direction;i++){
		salto_dir = salto_dir * Size[i];
	}


	for(i=0L;i<Size[direction];i++){
		Image[i*salto_dir+jump] = Line[i];
	}

}
/*--------------------------------------------------------------------------*/
extern int	SamplesToCoefficientsND
				(
					double	*Image,			/* in-place processing */
					long	*Size,			/* Vector Size of the Image */
					long	SplineDegree,	/* degree of the spline model */
					long	NDims			/* Number of dimension of the image*/
				){

	double	*Line;
	double	Pole[4];
	long	NbPoles;
	long	x, y,z,t;
	long coordinate[4];

	/* recover the poles from a lookup table */
	switch (SplineDegree) {
		case 0L:
			return 0;
			break;
		case 1L:
			return 0;
			break;
		case 2L:
			NbPoles = 1L;
			Pole[0] = sqrt(8.0) - 3.0;
			break;
		case 3L:
			NbPoles = 1L;
			Pole[0] = sqrt(3.0) - 2.0;
			break;
		case 4L:
			NbPoles = 2L;
			Pole[0] = sqrt(664.0 - sqrt(438976.0)) + sqrt(304.0) - 19.0;
			Pole[1] = sqrt(664.0 + sqrt(438976.0)) - sqrt(304.0) - 19.0;
			break;
		case 5L:
			NbPoles = 2L;
			Pole[0] = sqrt(135.0 / 2.0 - sqrt(17745.0 / 4.0)) + sqrt(105.0 / 4.0)
				- 13.0 / 2.0;
			Pole[1] = sqrt(135.0 / 2.0 + sqrt(17745.0 / 4.0)) - sqrt(105.0 / 4.0)
				- 13.0 / 2.0;
			break;
		case 6L:
			NbPoles = 3L;
			Pole[0] = -0.48829458930304475513011803888378906211227916123938;
			Pole[1] = -0.081679271076237512597937765737059080653379610398148;
			Pole[2] = -0.0014141518083258177510872439765585925278641690553467;
			break;
		case 7L:
			NbPoles = 3L;
			Pole[0] = -0.53528043079643816554240378168164607183392315234269;
			Pole[1] = -0.12255461519232669051527226435935734360548654942730;
			Pole[2] = -0.0091486948096082769285930216516478534156925639545994;
			break;
		case 8L:
			NbPoles = 4L;
			Pole[0] = -0.57468690924876543053013930412874542429066157804125;
			Pole[1] = -0.16303526929728093524055189686073705223476814550830;
			Pole[2] = -0.023632294694844850023403919296361320612665920854629;
			Pole[3] = -0.00015382131064169091173935253018402160762964054070043;
			break;
		case 9L:
			NbPoles = 4L;
			Pole[0] = -0.60799738916862577900772082395428976943963471853991;
			Pole[1] = -0.20175052019315323879606468505597043468089886575747;
			Pole[2] = -0.043222608540481752133321142979429688265852380231497;
			Pole[3] = -0.0021213069031808184203048965578486234220548560988624;
			break;
        case -3L: // oMom degree 3
            NbPoles = 1L;
			Pole[0] = -(1.0/8.0)*(13.0-sqrt(105.0));
            break;
		default:
			printf("Invalid spline degree\n");
			return(1);
	}


	switch(NDims){
	case 1:
        printf("Entro para convertirlos\n");
		ConvertToInterpolationCoefficients(Image, Size[0], Pole, NbPoles, DBL_EPSILON);
		break;
	case 2:
		// Convert x direction
		coordinate[0] = 0;
		for(y=0; y<Size[1]; y++){
			coordinate[1] = y;
			Line = GetMetaRow(Image,Size,coordinate,0,NDims);
			ConvertToInterpolationCoefficients(Line, Size[0], Pole, NbPoles, DBL_EPSILON);
			PutMetaRow(Image,Size,coordinate,0,NDims,Line);
			free(Line);
		}

		// Convert y direction
		coordinate[1] = 0;
		for(x=0; x<Size[0]; x++){
			coordinate[0] = x;
			Line = GetMetaRow(Image,Size,coordinate,1,NDims);
			ConvertToInterpolationCoefficients(Line, Size[1], Pole, NbPoles, DBL_EPSILON);
			PutMetaRow(Image,Size,coordinate,1,NDims,Line);
			free(Line);
		}
		break;
	case 3:
		// Convert x direction
		coordinate[0] = 0;
		for(y=0; y<Size[1]; y++){
			coordinate[1] = y;
			for(z=0;z<Size[2];z++){
				coordinate[2] = z;
				Line = GetMetaRow(Image,Size,coordinate,0,NDims);
				ConvertToInterpolationCoefficients(Line, Size[0], Pole, NbPoles, DBL_EPSILON);
				PutMetaRow(Image,Size,coordinate,0,NDims,Line);
				free(Line);
			}
		}
		// Convert y direction
		coordinate[1] = 0;
		for(x=0; x<Size[0]; x++){
			coordinate[0] = x;
			for(z=0;z<Size[2];z++){
				coordinate[2] = z;
				Line = GetMetaRow(Image,Size,coordinate,1,NDims);
				ConvertToInterpolationCoefficients(Line, Size[1], Pole, NbPoles, DBL_EPSILON);
				PutMetaRow(Image,Size,coordinate,1,NDims,Line);
				free(Line);
			}
		}
		// Convert z direction
		coordinate[2] = 0;
		for(x=0; x<Size[0]; x++){
			coordinate[0] = x;
			for(y=0;y<Size[1];y++){
				coordinate[1] = y;
				Line = GetMetaRow(Image,Size,coordinate,2,NDims);
				ConvertToInterpolationCoefficients(Line, Size[2], Pole, NbPoles, DBL_EPSILON);
				PutMetaRow(Image,Size,coordinate,2,NDims,Line);
				free(Line);
			}
		}
		break;
	case 4:
		// Convert x direction
		coordinate[0] = 0;
		for(y=0; y<Size[1]; y++){
			coordinate[1] = y;
			for(z=0;z<Size[2];z++){
				coordinate[2] = z;
				for(t=0;t<Size[3];t++){
					coordinate[3] = t;
					Line = GetMetaRow(Image,Size,coordinate,0,NDims);
					ConvertToInterpolationCoefficients(Line, Size[0], Pole, NbPoles, DBL_EPSILON);
					PutMetaRow(Image,Size,coordinate,0,NDims,Line);
					free(Line);
				}
			}
		}
		// Convert y direction
		coordinate[1] = 0;
		for(x=0; x<Size[0]; x++){
			coordinate[0] = x;
			for(z=0;z<Size[2];z++){
				coordinate[2] = z;
				for(t=0;t<Size[3];t++){
					coordinate[3] = t;
					Line = GetMetaRow(Image,Size,coordinate,1,NDims);
					ConvertToInterpolationCoefficients(Line, Size[1], Pole, NbPoles, DBL_EPSILON);
					PutMetaRow(Image,Size,coordinate,1,NDims,Line);
					free(Line);
				}
			}
		}
		// Convert z direction
		coordinate[2] = 0;
		for(x=0; x<Size[0]; x++){
			coordinate[0] = x;
			for(y=0;y<Size[1];y++){
				coordinate[1] = y;
				for(t=0;t<Size[3];t++){
					coordinate[3] = t;
					Line = GetMetaRow(Image,Size,coordinate,2,NDims);
					ConvertToInterpolationCoefficients(Line, Size[2], Pole, NbPoles, DBL_EPSILON);
					PutMetaRow(Image,Size,coordinate,2,NDims,Line);
					free(Line);
				}
			}
		}
		// Convert t direction
		coordinate[3] = 0;
		for(x=0; x<Size[0]; x++){
			coordinate[0] = x;
			for(y=0;y<Size[1];y++){
				coordinate[1] = y;
				for(z=0;z<Size[2];z++){
					coordinate[2] = z;
					Line = GetMetaRow(Image,Size,coordinate,3,NDims);
					ConvertToInterpolationCoefficients(Line, Size[3], Pole, NbPoles, DBL_EPSILON);
					PutMetaRow(Image,Size,coordinate,3,NDims,Line);
					free(Line);
				}
			}
		}
		break;
	}


return 0;


}



/*****************************************************************************
 *	Definition of extern procedures
 ****************************************************************************/
/*--------------------------------------------------------------------------*/
extern int		SamplesToCoefficients
				(
					float	*Image,		/* in-place processing */
					long	Width,		/* width of the image */
					long	Height,		/* height of the image */
					long	SplineDegree/* degree of the spline model */
				)

{ /* begin SamplesToCoefficients */

	double	*Line;
	double	Pole[4];
	long	NbPoles;
	long	x, y;

	/* recover the poles from a lookup table */
	switch (SplineDegree) {
		case 2L:
			NbPoles = 1L;
			Pole[0] = sqrt(8.0) - 3.0;
			break;
		case 3L:
			NbPoles = 1L;
			Pole[0] = sqrt(3.0) - 2.0;
			break;
		case 4L:
			NbPoles = 2L;
			Pole[0] = sqrt(664.0 - sqrt(438976.0)) + sqrt(304.0) - 19.0;
			Pole[1] = sqrt(664.0 + sqrt(438976.0)) - sqrt(304.0) - 19.0;
			break;
		case 5L:
			NbPoles = 2L;
			Pole[0] = sqrt(135.0 / 2.0 - sqrt(17745.0 / 4.0)) + sqrt(105.0 / 4.0)
				- 13.0 / 2.0;
			Pole[1] = sqrt(135.0 / 2.0 + sqrt(17745.0 / 4.0)) - sqrt(105.0 / 4.0)
				- 13.0 / 2.0;
			break;
		case 6L:
			NbPoles = 3L;
			Pole[0] = -0.48829458930304475513011803888378906211227916123938;
			Pole[1] = -0.081679271076237512597937765737059080653379610398148;
			Pole[2] = -0.0014141518083258177510872439765585925278641690553467;
			break;
		case 7L:
			NbPoles = 3L;
			Pole[0] = -0.53528043079643816554240378168164607183392315234269;
			Pole[1] = -0.12255461519232669051527226435935734360548654942730;
			Pole[2] = -0.0091486948096082769285930216516478534156925639545994;
			break;
		case 8L:
			NbPoles = 4L;
			Pole[0] = -0.57468690924876543053013930412874542429066157804125;
			Pole[1] = -0.16303526929728093524055189686073705223476814550830;
			Pole[2] = -0.023632294694844850023403919296361320612665920854629;
			Pole[3] = -0.00015382131064169091173935253018402160762964054070043;
			break;
		case 9L:
			NbPoles = 4L;
			Pole[0] = -0.60799738916862577900772082395428976943963471853991;
			Pole[1] = -0.20175052019315323879606468505597043468089886575747;
			Pole[2] = -0.043222608540481752133321142979429688265852380231497;
			Pole[3] = -0.0021213069031808184203048965578486234220548560988624;
			break;
		default:
			printf("Invalid spline degree\n");
			return(1);
	}

	/* convert the image samples into interpolation coefficients */
	/* in-place separable process, along x */
	Line = (double *)malloc((size_t)(Width * (long)sizeof(double)));
	if (Line == (double *)NULL) {
		printf("Row allocation failed\n");
		return(1);
	}
	for (y = 0L; y < Height; y++) {
		GetRow(Image, y, Line, Width);
		ConvertToInterpolationCoefficients(Line, Width, Pole, NbPoles, DBL_EPSILON);
		PutRow(Image, y, Line, Width);
	}
	free(Line);

	/* in-place separable process, along y */
	Line = (double *)malloc((size_t)(Height * (long)sizeof(double)));
	if (Line == (double *)NULL) {
		printf("Column allocation failed\n");
		return(1);
	}
	for (x = 0L; x < Width; x++) {
		GetColumn(Image, Width, x, Line, Height);
		ConvertToInterpolationCoefficients(Line, Height, Pole, NbPoles, DBL_EPSILON);
		PutColumn(Image, Width, x, Line, Height);
	}
	free(Line);

	return(0);
} /* end SamplesToCoefficients */





