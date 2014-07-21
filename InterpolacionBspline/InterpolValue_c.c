/****************************************************************************
 *
 * Date: February 20, 2013
 *
 * E. Marinetto
 *
 * This function allow to use bspline_interpol software in Matlab software.
 *
 *
 *
 * E. Marinetto
 * Laboratorio de Imagen Médica
 * WEB del LIM
 * Unidad de Medicina y Cirugía Experimental.
 * Fundación para la Investigación Biomédica del
 * Hospital Gregorio Marañón
 * C/ Doctor Esquerdo, 46
 * 28007  Madrid
 * Teléfono: +34 91 4265017
 * Fax: +34 91 4265108
 *
****************************************************************************/

/*****************************************************************************
 *	Include
 ****************************************************************************/
#include	"InterpolValue_c.h"

/*****************************************************************************
 *	This file implements the method for mexfile
 *
 *	It convert a vector of data to the coefficients values.
 *	Input is passed as double
 ****************************************************************************/



/*************************************************************************/


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){



	/* This function implements a interface to matlab in order to use :
	 *
	 * extern int		SamplesToCoefficients
				(
					double	*Line,		/* in-place processing *
					int	L,	/* dimension of the image *
					int	SplineDegree	/* degree of the spline model *
				);
	 */

	/* Inputs */
    #define COEFFS		prhs[0]
    #define COORDINATE	prhs[1]
    #define SPLINEDEG	prhs[2]

	/* Outputs */
    #define INTERPOLATEDVALUE  	plhs[0]      
            
	/* Duplicate the input to get the output */

    double* Bcoeff;
    double* coordinates;
    long    NDims;
    long    SplineDegree;
    double* valueInterpolated;
    double  coordinate[4];
    long    Size[4];
    long*   size_aux;
    int     iterator;
    double  valor;
    
    
    if(!mxIsDouble(COEFFS)){
        mexErrMsgTxt("La matriz de coeficientes tiene que ser tipo double");
    }  
    
    Bcoeff = mxGetPr(COEFFS);
    coordinates = mxGetPr(COORDINATE);
    NDims = (long) mxGetNumberOfElements(COORDINATE);
    SplineDegree = (long) mxGetScalar(SPLINEDEG);
    size_aux = (long*)mxGetDimensions(COEFFS);
    
    INTERPOLATEDVALUE = mxCreateDoubleMatrix(1,1,mxREAL);
    valueInterpolated = (double *) mxGetPr(INTERPOLATEDVALUE);
    valueInterpolated[0] = 0;

	for(iterator=0; iterator<4; iterator++){
		 coordinate[iterator] = 0;
			 Size[iterator] = 1;
	}

    
    
	for(iterator=0;iterator<mxGetNumberOfDimensions(COEFFS);iterator++){
		 coordinate[iterator] = coordinates[iterator];
		 Size[iterator] = (long) size_aux[iterator];
         
	}
    

    InterpolatedValueND_matlab(Bcoeff,coordinate,Size,NDims,SplineDegree,&valor);
    valueInterpolated[0] = valor;
    return;
}


