/****************************************************************************
 *
 * Date: February 20, 2013
 *
 * E. Marinetto
 *
 * This function allow to use bspline_coeff software in Matlab software.
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
#include	"ConvertToInterpolationCoefficients_c.h"

/*****************************************************************************
 *	This file implements the method for mexfile
 *
 *	It convert a vector of data to the coefficients values.
 *	Input is passed as double up to 4 dimensions
 ****************************************************************************/

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
    #define SAMPLES		prhs[0]
    #define SPLINEDEG	prhs[1]

	/* Outputs */
    #define COEFFS  	plhs[0]
            
	mxArray* output;
    double* samples_copy;
	int SplineDegree;
	int error;
    long NDims;
    long* Size;
    
    
    if(!mxIsDouble(SAMPLES)){
        mexErrMsgTxt("La matriz de entrada tiene que ser tipo double");
    }
    
    output = mxDuplicateArray(SAMPLES);
    samples_copy = (double*)mxGetPr(output);
	SplineDegree = (int)mxGetScalar(SPLINEDEG);
    NDims = mxGetNumberOfDimensions(SAMPLES);
    Size = (long*) mxGetDimensions(SAMPLES);
    
    
    error = SamplesToCoefficientsND(samples_copy,Size,SplineDegree,NDims);

    COEFFS = output;

}
