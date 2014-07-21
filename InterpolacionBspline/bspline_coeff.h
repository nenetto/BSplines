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
static void		ConvertToInterpolationCoefficients
				(
					double	c[],		/* input samples --> output coefficients */
					long	DataLength,	/* number of samples or coefficients */
					double	z[],		/* poles */
					long	NbPoles,	/* number of poles */
					double	Tolerance	/* admissible relative error */
				);

/*--------------------------------------------------------------------------*/
static void		GetColumn
				(
					float	*Image,		/* input image array */
					long	Width,		/* width of the image */
					long	x,			/* x coordinate of the selected line */
					double	Line[],		/* output linear array */
					long	Height		/* length of the line and height of the image */
				);

/*--------------------------------------------------------------------------*/
static void		GetRow
				(
					float	*Image,		/* input image array */
					long	y,			/* y coordinate of the selected line */
					double	Line[],		/* output linear array */
					long	Width		/* length of the line and width of the image */
				);

/*--------------------------------------------------------------------------*/
static double	InitialCausalCoefficient
				(
					double	c[],		/* coefficients */
					long	DataLength,	/* number of coefficients */
					double	z,			/* actual pole */
					double	Tolerance	/* admissible relative error */
				);

/*--------------------------------------------------------------------------*/
static double	InitialAntiCausalCoefficient
				(
					double	c[],		/* coefficients */
					long	DataLength,	/* number of samples or coefficients */
					double	z			/* actual pole */
				);

/*--------------------------------------------------------------------------*/
static void		PutColumn
				(
					float	*Image,		/* output image array */
					long	Width,		/* width of the image */
					long	x,			/* x coordinate of the selected line */
					double	Line[],		/* input linear array */
					long	Height		/* length of the line and height of the image */
				);

/*--------------------------------------------------------------------------*/
static void		PutRow
				(
					float	*Image,		/* output image array */
					long	y,			/* y coordinate of the selected line */
					double	Line[],		/* input linear array */
					long	Width		/* length of the line and width of the image */
				);
/*--------------------------------------------------------------------------*/
/*****************************************************************************
 *	N-D Methods by E. Marinetto
 ****************************************************************************/

static double*	GetMetaRow
				(
					double	*Image,			/* input image array */
					long	*Size,			/* array with the dimensions of the image */
					long	*coordinate,	/* coordinate of the selected line */
					long	direction,		/* direction of the line */
					long	ImageDim		/* number of dimensions of the image */
				);

static void 	PutMetaRow
				(
					double	*Image,			/* input image array */
					long	*Size,			/* array with the dimensions of the image */
					long	*coordinate,	/* coordinate of the selected line */
					long	direction,		/* direction of the row */
					long	ImageDim,		/* number of dimensions of the image */
					double  Line[]			/* line to put in*/
				);

extern int		SamplesToCoefficientsND
				(
					double	*Image,			/* in-place processing */
					long	*Size,			/* Vector Size of the Image */
					long	SplineDegree,	/* degree of the spline model */
					long	NDims			/* number of dimension of the image*/
				);






