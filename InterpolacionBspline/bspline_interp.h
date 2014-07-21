/*****************************************************************************
 * *
 *	Date: February 20, 2013
 *
 *	Some modifications over previous version were made to accomplish a
 *	N-Dimensional interpolation. Currently, this software produce Bspline
 *	interpolation and oMom3 interpolation  up to 4 dimensions.
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
 *	System includes
 ****************************************************************************/
#include	<float.h>
#include	<math.h>
#include	<stddef.h>
#include	<stdio.h>
#include	<stdlib.h>
#define printf mexPrintf
//#include    "mex.h"

/*****************************************************************************
 *	Declaration of static procedures
 ****************************************************************************/
/*--------------------------------------------------------------------------*/

extern double	InterpolatedValueND
				(
					double	*Bcoeff,		/* input B-spline array of coefficients */
					double	*coordinate,	/*  coordinate where to interpolate */
                    long    *Size,      	/* Size of the image */
					long	NDims,			/* y coordinate where to interpolate */
					long	SplineDegree	/* degree of the spline model */
				);

extern void	InterpolatedValueND_matlab
				(
					double	*Bcoeff,		/* input B-spline array of coefficients */
					double	*coordinate,	/*  coordinate where to interpolate */
                    long    *Size,      	/* Size of the image */
					long	NDims,			/* y coordinate where to interpolate */
					long	SplineDegree,	/* degree of the spline model */
                    double  *interpolated
				);

extern double	oMom3(
					double x				/*input for the oMom3 kernel */
				);
extern double lin(
					double x				/*input for the lineal kernel */
				);
