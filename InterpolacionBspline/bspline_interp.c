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
 *	Includes
 ****************************************************************************/

#include "bspline_interp.h"


/*****************************************************************************
 *	Definition of static procedures
 ****************************************************************************/
/*--------------------------------------------------------------------------*/

extern void	InterpolatedValueND_matlab
				(
					double	*Bcoeff,		/* input B-spline array of coefficients */
					double	*coordinate,	/*  coordinate where to interpolate */
                    long    *Size,      	/* Size of the image */
					long	NDims,			/* y coordinate where to interpolate */
					long	SplineDegree,	/* degree of the spline model */
                    double  *interpolated
				){
 
    double valor = InterpolatedValueND(Bcoeff,coordinate,Size,NDims,SplineDegree);
    interpolated[0] = valor;
    
    
}


extern double	InterpolatedValueND
				(
					double	*Bcoeff,	/* input B-spline array of coefficients */
					double	*coordinate,	/*  coordinate where to interpolate */
                    long    *Size,      /* Size of the image */
					long	NDims,		/* y coordinate where to interpolate */
					long	SplineDegree/* degree of the spline model */
				){


    
   
	double	xWeight[10], yWeight[10], zWeight[10], tWeight[10];
	double	interpolated;
	double	w, w2, w4, t, t0, t1,w3;
	long	xIndex[10], yIndex[10], zIndex[10], tIndex[10];
    long    size_X = Size[0];
    long    size_Y = Size[1];
    long    size_Z = Size[2];
    long    size_T = Size[3];
	long	size_X2 = 2L * Size[0] - 2L;
    long    size_Y2 = 2L * Size[1] - 2L;
    long    size_Z2 = 2L * Size[2] - 2L;
    long    size_T2 = 2L * Size[3] - 2L;
	long	i, j, i2,j2,k,n;
    
    //printf("Size1 = %ld %ld %ld %ld\n",Size[0],Size[1],Size[2],Size[3]);
    //printf("SplineDegree = %ld\n",SplineDegree);
 	/* compute the interpolation indexes */
 	if (abs(SplineDegree) & 1L) {
			i = (long)floor(coordinate[0]) - abs(SplineDegree) / 2L;
			j = (long)floor(coordinate[1]) - abs(SplineDegree) / 2L;
			i2 = (long)floor(coordinate[2]) - abs(SplineDegree) / 2L;
			j2 = (long)floor(coordinate[3]) - abs(SplineDegree) / 2L;
        //printf("ENTRADA A\n");
 		for (k = 0L; k <= abs(SplineDegree); k++) {
			xIndex[k] = i++;
			yIndex[k] = j++;
			zIndex[k] = i2++;
			tIndex[k] = j2++;
            //printf("ENTRADA B\n");
 		}
 	}
 	else {
        //printf("ENTRADA C\n");
		i = (long)floor(coordinate[0] + 0.5) - abs(SplineDegree) / 2L;
		j = (long)floor(coordinate[1] + 0.5) - abs(SplineDegree) / 2L;
		i2 = (long)floor(coordinate[2] + 0.5) - abs(SplineDegree) / 2L;
		j2 = (long)floor(coordinate[3] + 0.5) - abs(SplineDegree) / 2L;
 		for (k = 0L; k <= abs(SplineDegree); k++) {
            
			xIndex[k] = i++;
			yIndex[k] = j++;
			zIndex[k] = i2++;
			tIndex[k] = j2++;
 		}
 	}

    
    
    
    //printf("Size2 = %ld %ld %ld %ld\n",Size[0],Size[1],Size[2],Size[3]);
    
	/* compute the interpolation weights */
 	switch (SplineDegree) {
		case 0L:
			/* x */
			xWeight[0] = 1;
			/* y */
			yWeight[0] = 1;
			/* z */
			zWeight[0] = 1;
			/* t */
			tWeight[0] = 1;
			break;
		case 1L:
			/* x */
            xWeight[0] = lin(coordinate[0] - (double)xIndex[0]);
			xWeight[1] = 1.0-xWeight[0];
			/* y */
			yWeight[0] = lin(coordinate[1] - (double)yIndex[0]);
			yWeight[1] = 1.0-yWeight[0];
			/* z */
			zWeight[0] = lin(coordinate[2] - (double)zIndex[0]);
			zWeight[1] = 1.0-zWeight[0];
			/* t */
			tWeight[0] = lin(coordinate[3] - (double)tIndex[0]);
			tWeight[1] = 1.0-tWeight[0];
           //printf("Size3 = %ld %ld %ld %ld\n",Size[0],Size[1],Size[2],Size[3]);
            break;
 		case 2L:
 			/* x */
 			w = coordinate[0] - (double)xIndex[1];
 			xWeight[1] = 3.0 / 4.0 - w * w;
 			xWeight[2] = (1.0 / 2.0) * (w - xWeight[1] + 1.0);
 			xWeight[0] = 1.0 - xWeight[1] - xWeight[2];
 			/* y */
 			w = coordinate[1] - (double)yIndex[1];
 			yWeight[1] = 3.0 / 4.0 - w * w;
 			yWeight[2] = (1.0 / 2.0) * (w - yWeight[1] + 1.0);
 			yWeight[0] = 1.0 - yWeight[1] - yWeight[2];
             /* z */
 			w = coordinate[2] - (double)zIndex[1];
 			zWeight[1] = 3.0 / 4.0 - w * w;
 			zWeight[2] = (1.0 / 2.0) * (w - zWeight[1] + 1.0);
 			zWeight[0] = 1.0 - zWeight[1] - zWeight[2];
             /* t */
 			w = coordinate[3] - (double)tIndex[1];
 			tWeight[1] = 3.0 / 4.0 - w * w;
 			tWeight[2] = (1.0 / 2.0) * (w - tWeight[1] + 1.0);
 			tWeight[0] = 1.0 - tWeight[1] - tWeight[2];
 			break;
 		case 3L:
 			/* x */
 			w = coordinate[0] - (double)xIndex[1];
 			xWeight[3] = (1.0 / 6.0) * w * w * w;
 			xWeight[0] = (1.0 / 6.0) + (1.0 / 2.0) * w * (w - 1.0) - xWeight[3];
 			xWeight[2] = w + xWeight[0] - 2.0 * xWeight[3];
 			xWeight[1] = 1.0 - xWeight[0] - xWeight[2] - xWeight[3];
 			/* y */
 			w = coordinate[1] - (double)yIndex[1];
 			yWeight[3] = (1.0 / 6.0) * w * w * w;
 			yWeight[0] = (1.0 / 6.0) + (1.0 / 2.0) * w * (w - 1.0) - yWeight[3];
 			yWeight[2] = w + yWeight[0] - 2.0 * yWeight[3];
 			yWeight[1] = 1.0 - yWeight[0] - yWeight[2] - yWeight[3];
             /* z */
 			w = coordinate[2] - (double)zIndex[1];
 			zWeight[3] = (1.0 / 6.0) * w * w * w;
 			zWeight[0] = (1.0 / 6.0) + (1.0 / 2.0) * w * (w - 1.0) - zWeight[3];
 			zWeight[2] = w + zWeight[0] - 2.0 * zWeight[3];
 			zWeight[1] = 1.0 - zWeight[0] - zWeight[2] - zWeight[3];
             /* t */
 			w = coordinate[3] - (double)tIndex[1];
 			tWeight[3] = (1.0 / 6.0) * w * w * w;
 			tWeight[0] = (1.0 / 6.0) + (1.0 / 2.0) * w * (w - 1.0) - tWeight[3];
 			tWeight[2] = w + tWeight[0] - 2.0 * tWeight[3];
 			tWeight[1] = 1.0 - tWeight[0] - tWeight[2] - tWeight[3];
 			break;
 		case 4L:
 			/* x */
 			w = coordinate[0] - (double)xIndex[2];
 			w2 = w * w;
 			t = (1.0 / 6.0) * w2;
 			xWeight[0] = 1.0 / 2.0 - w;
 			xWeight[0] *= xWeight[0];
 			xWeight[0] *= (1.0 / 24.0) * xWeight[0];
 			t0 = w * (t - 11.0 / 24.0);
 			t1 = 19.0 / 96.0 + w2 * (1.0 / 4.0 - t);
 			xWeight[1] = t1 + t0;
 			xWeight[3] = t1 - t0;
 			xWeight[4] = xWeight[0] + t0 + (1.0 / 2.0) * w;
 			xWeight[2] = 1.0 - xWeight[0] - xWeight[1] - xWeight[3] - xWeight[4];
 			/* y */
 			w = coordinate[1] - (double)yIndex[2];
 			w2 = w * w;
 			t = (1.0 / 6.0) * w2;
 			yWeight[0] = 1.0 / 2.0 - w;
 			yWeight[0] *= yWeight[0];
 			yWeight[0] *= (1.0 / 24.0) * yWeight[0];
 			t0 = w * (t - 11.0 / 24.0);
 			t1 = 19.0 / 96.0 + w2 * (1.0 / 4.0 - t);
 			yWeight[1] = t1 + t0;
 			yWeight[3] = t1 - t0;
 			yWeight[4] = yWeight[0] + t0 + (1.0 / 2.0) * w;
 			yWeight[2] = 1.0 - yWeight[0] - yWeight[1] - yWeight[3] - yWeight[4];
             /* z */
 			w = coordinate[2] - (double)zIndex[2];
 			w2 = w * w;
 			t = (1.0 / 6.0) * w2;
 			zWeight[0] = 1.0 / 2.0 - w;
 			zWeight[0] *= zWeight[0];
 			zWeight[0] *= (1.0 / 24.0) * zWeight[0];
 			t0 = w * (t - 11.0 / 24.0);
 			t1 = 19.0 / 96.0 + w2 * (1.0 / 4.0 - t);
 			zWeight[1] = t1 + t0;
 			zWeight[3] = t1 - t0;
 			zWeight[4] = zWeight[0] + t0 + (1.0 / 2.0) * w;
 			zWeight[2] = 1.0 - zWeight[0] - zWeight[1] - zWeight[3] - zWeight[4];
             /* t */
 			w = coordinate[3] - (double)tIndex[2];
 			w2 = w * w;
 			t = (1.0 / 6.0) * w2;
 			tWeight[0] = 1.0 / 2.0 - w;
 			tWeight[0] *= tWeight[0];
 			tWeight[0] *= (1.0 / 24.0) * tWeight[0];
 			t0 = w * (t - 11.0 / 24.0);
 			t1 = 19.0 / 96.0 + w2 * (1.0 / 4.0 - t);
 			tWeight[1] = t1 + t0;
 			tWeight[3] = t1 - t0;
 			tWeight[4] = tWeight[0] + t0 + (1.0 / 2.0) * w;
 			tWeight[2] = 1.0 - tWeight[0] - tWeight[1] - tWeight[3] - tWeight[4];

 			break;
 		case 5L:
 			/* x */
 			w = coordinate[0] - (double)xIndex[2];
 			w2 = w * w;
 			xWeight[5] = (1.0 / 120.0) * w * w2 * w2;
 			w2 -= w;
 			w4 = w2 * w2;
 			w -= 1.0 / 2.0;
 			t = w2 * (w2 - 3.0);
 			xWeight[0] = (1.0 / 24.0) * (1.0 / 5.0 + w2 + w4) - xWeight[5];
 			t0 = (1.0 / 24.0) * (w2 * (w2 - 5.0) + 46.0 / 5.0);
 			t1 = (-1.0 / 12.0) * w * (t + 4.0);
 			xWeight[2] = t0 + t1;
 			xWeight[3] = t0 - t1;
 			t0 = (1.0 / 16.0) * (9.0 / 5.0 - t);
 			t1 = (1.0 / 24.0) * w * (w4 - w2 - 5.0);
 			xWeight[1] = t0 + t1;
 			xWeight[4] = t0 - t1;
 			/* y */
 			w = coordinate[1] - (double)yIndex[2];
 			w2 = w * w;
 			yWeight[5] = (1.0 / 120.0) * w * w2 * w2;
 			w2 -= w;
 			w4 = w2 * w2;
 			w -= 1.0 / 2.0;
 			t = w2 * (w2 - 3.0);
 			yWeight[0] = (1.0 / 24.0) * (1.0 / 5.0 + w2 + w4) - yWeight[5];
 			t0 = (1.0 / 24.0) * (w2 * (w2 - 5.0) + 46.0 / 5.0);
 			t1 = (-1.0 / 12.0) * w * (t + 4.0);
 			yWeight[2] = t0 + t1;
 			yWeight[3] = t0 - t1;
 			t0 = (1.0 / 16.0) * (9.0 / 5.0 - t);
 			t1 = (1.0 / 24.0) * w * (w4 - w2 - 5.0);
 			yWeight[1] = t0 + t1;
 			yWeight[4] = t0 - t1;
             /* z */
 			w = coordinate[2] - (double)zIndex[2];
 			w2 = w * w;
 			zWeight[5] = (1.0 / 120.0) * w * w2 * w2;
 			w2 -= w;
 			w4 = w2 * w2;
 			w -= 1.0 / 2.0;
 			t = w2 * (w2 - 3.0);
 			zWeight[0] = (1.0 / 24.0) * (1.0 / 5.0 + w2 + w4) - zWeight[5];
 			t0 = (1.0 / 24.0) * (w2 * (w2 - 5.0) + 46.0 / 5.0);
 			t1 = (-1.0 / 12.0) * w * (t + 4.0);
 			zWeight[2] = t0 + t1;
 			zWeight[3] = t0 - t1;
 			t0 = (1.0 / 16.0) * (9.0 / 5.0 - t);
 			t1 = (1.0 / 24.0) * w * (w4 - w2 - 5.0);
 			zWeight[1] = t0 + t1;
 			zWeight[4] = t0 - t1;
             /* t */
 			w = coordinate[0] - (double)tIndex[2];
 			w2 = w * w;
 			tWeight[5] = (1.0 / 120.0) * w * w2 * w2;
 			w2 -= w;
 			w4 = w2 * w2;
 			w -= 1.0 / 2.0;
 			t = w2 * (w2 - 3.0);
 			tWeight[0] = (1.0 / 24.0) * (1.0 / 5.0 + w2 + w4) - tWeight[5];
 			t0 = (1.0 / 24.0) * (w2 * (w2 - 5.0) + 46.0 / 5.0);
 			t1 = (-1.0 / 12.0) * w * (t + 4.0);
 			tWeight[2] = t0 + t1;
 			tWeight[3] = t0 - t1;
 			t0 = (1.0 / 16.0) * (9.0 / 5.0 - t);
 			t1 = (1.0 / 24.0) * w * (w4 - w2 - 5.0);
 			tWeight[1] = t0 + t1;
 			tWeight[4] = t0 - t1;
 			break;
 		case 6L:
 			/* x */
 			w = coordinate[0] - (double)xIndex[3];
 			xWeight[0] = 1.0 / 2.0 - w;
 			xWeight[0] *= xWeight[0] * xWeight[0];
 			xWeight[0] *= xWeight[0] / 720.0;
 			xWeight[1] = (361.0 / 192.0 - w * (59.0 / 8.0 + w
 				* (-185.0 / 16.0 + w * (25.0 / 3.0 + w * (-5.0 / 2.0 + w)
 				* (1.0 / 2.0 + w))))) / 120.0;
 			xWeight[2] = (10543.0 / 960.0 + w * (-289.0 / 16.0 + w
 				* (79.0 / 16.0 + w * (43.0 / 6.0 + w * (-17.0 / 4.0 + w
 				* (-1.0 + w)))))) / 48.0;
 			w2 = w * w;
 			xWeight[3] = (5887.0 / 320.0 - w2 * (231.0 / 16.0 - w2
 				* (21.0 / 4.0 - w2))) / 36.0;
 			xWeight[4] = (10543.0 / 960.0 + w * (289.0 / 16.0 + w
 				* (79.0 / 16.0 + w * (-43.0 / 6.0 + w * (-17.0 / 4.0 + w
 				* (1.0 + w)))))) / 48.0;
 			xWeight[6] = 1.0 / 2.0 + w;
 			xWeight[6] *= xWeight[6] * xWeight[6];
 			xWeight[6] *= xWeight[6] / 720.0;
 			xWeight[5] = 1.0 - xWeight[0] - xWeight[1] - xWeight[2] - xWeight[3]
 				- xWeight[4] - xWeight[6];
 			/* y */
 			w = coordinate[1] - (double)yIndex[3];
 			yWeight[0] = 1.0 / 2.0 - w;
 			yWeight[0] *= yWeight[0] * yWeight[0];
 			yWeight[0] *= yWeight[0] / 720.0;
 			yWeight[1] = (361.0 / 192.0 - w * (59.0 / 8.0 + w
 				* (-185.0 / 16.0 + w * (25.0 / 3.0 + w * (-5.0 / 2.0 + w)
 				* (1.0 / 2.0 + w))))) / 120.0;
 			yWeight[2] = (10543.0 / 960.0 + w * (-289.0 / 16.0 + w
 				* (79.0 / 16.0 + w * (43.0 / 6.0 + w * (-17.0 / 4.0 + w
 				* (-1.0 + w)))))) / 48.0;
 			w2 = w * w;
 			yWeight[3] = (5887.0 / 320.0 - w2 * (231.0 / 16.0 - w2
 				* (21.0 / 4.0 - w2))) / 36.0;
 			yWeight[4] = (10543.0 / 960.0 + w * (289.0 / 16.0 + w
 				* (79.0 / 16.0 + w * (-43.0 / 6.0 + w * (-17.0 / 4.0 + w
 				* (1.0 + w)))))) / 48.0;
 			yWeight[6] = 1.0 / 2.0 + w;
 			yWeight[6] *= yWeight[6] * yWeight[6];
 			yWeight[6] *= yWeight[6] / 720.0;
 			yWeight[5] = 1.0 - yWeight[0] - yWeight[1] - yWeight[2] - yWeight[3]
 				- yWeight[4] - yWeight[6];
             /* z */
 			w = coordinate[2] - (double)zIndex[3];
 			zWeight[0] = 1.0 / 2.0 - w;
 			zWeight[0] *= zWeight[0] * zWeight[0];
 			zWeight[0] *= zWeight[0] / 720.0;
 			zWeight[1] = (361.0 / 192.0 - w * (59.0 / 8.0 + w
 				* (-185.0 / 16.0 + w * (25.0 / 3.0 + w * (-5.0 / 2.0 + w)
 				* (1.0 / 2.0 + w))))) / 120.0;
 			zWeight[2] = (10543.0 / 960.0 + w * (-289.0 / 16.0 + w
 				* (79.0 / 16.0 + w * (43.0 / 6.0 + w * (-17.0 / 4.0 + w
 				* (-1.0 + w)))))) / 48.0;
 			w2 = w * w;
 			zWeight[3] = (5887.0 / 320.0 - w2 * (231.0 / 16.0 - w2
 				* (21.0 / 4.0 - w2))) / 36.0;
 			zWeight[4] = (10543.0 / 960.0 + w * (289.0 / 16.0 + w
 				* (79.0 / 16.0 + w * (-43.0 / 6.0 + w * (-17.0 / 4.0 + w
 				* (1.0 + w)))))) / 48.0;
 			zWeight[6] = 1.0 / 2.0 + w;
 			zWeight[6] *= zWeight[6] * zWeight[6];
 			zWeight[6] *= zWeight[6] / 720.0;
 			zWeight[5] = 1.0 - zWeight[0] - zWeight[1] - zWeight[2] - zWeight[3]
 				- zWeight[4] - zWeight[6];
             /* t */
 			w = coordinate[3] - (double)tIndex[3];
 			tWeight[0] = 1.0 / 2.0 - w;
 			tWeight[0] *= tWeight[0] * tWeight[0];
 			tWeight[0] *= tWeight[0] / 720.0;
 			tWeight[1] = (361.0 / 192.0 - w * (59.0 / 8.0 + w
 				* (-185.0 / 16.0 + w * (25.0 / 3.0 + w * (-5.0 / 2.0 + w)
 				* (1.0 / 2.0 + w))))) / 120.0;
 			tWeight[2] = (10543.0 / 960.0 + w * (-289.0 / 16.0 + w
 				* (79.0 / 16.0 + w * (43.0 / 6.0 + w * (-17.0 / 4.0 + w
 				* (-1.0 + w)))))) / 48.0;
 			w2 = w * w;
 			tWeight[3] = (5887.0 / 320.0 - w2 * (231.0 / 16.0 - w2
 				* (21.0 / 4.0 - w2))) / 36.0;
 			tWeight[4] = (10543.0 / 960.0 + w * (289.0 / 16.0 + w
 				* (79.0 / 16.0 + w * (-43.0 / 6.0 + w * (-17.0 / 4.0 + w
 				* (1.0 + w)))))) / 48.0;
 			tWeight[6] = 1.0 / 2.0 + w;
 			tWeight[6] *= tWeight[6] * tWeight[6];
 			tWeight[6] *= tWeight[6] / 720.0;
 			tWeight[5] = 1.0 - tWeight[0] - tWeight[1] - tWeight[2] - tWeight[3]
 				- tWeight[4] - tWeight[6];
 			break;
 		case 7L:
 			/* x */
 			w = coordinate[0] - (double)xIndex[3];
 			xWeight[0] = 1.0 - w;
 			xWeight[0] *= xWeight[0];
 			xWeight[0] *= xWeight[0] * xWeight[0];
 			xWeight[0] *= (1.0 - w) / 5040.0;
 			w2 = w * w;
 			xWeight[1] = (120.0 / 7.0 + w * (-56.0 + w * (72.0 + w
 				* (-40.0 + w2 * (12.0 + w * (-6.0 + w)))))) / 720.0;
 			xWeight[2] = (397.0 / 7.0 - w * (245.0 / 3.0 + w * (-15.0 + w
 				* (-95.0 / 3.0 + w * (15.0 + w * (5.0 + w
 				* (-5.0 + w))))))) / 240.0;
 			xWeight[3] = (2416.0 / 35.0 + w2 * (-48.0 + w2 * (16.0 + w2
 				* (-4.0 + w)))) / 144.0;
 			xWeight[4] = (1191.0 / 35.0 - w * (-49.0 + w * (-9.0 + w
 				* (19.0 + w * (-3.0 + w) * (-3.0 + w2))))) / 144.0;
 			xWeight[5] = (40.0 / 7.0 + w * (56.0 / 3.0 + w * (24.0 + w
 				* (40.0 / 3.0 + w2 * (-4.0 + w * (-2.0 + w)))))) / 240.0;
 			xWeight[7] = w2;
 			xWeight[7] *= xWeight[7] * xWeight[7];
 			xWeight[7] *= w / 5040.0;
 			xWeight[6] = 1.0 - xWeight[0] - xWeight[1] - xWeight[2] - xWeight[3]
 				- xWeight[4] - xWeight[5] - xWeight[7];
 			/* y */
 			w = coordinate[1] - (double)yIndex[3];
 			yWeight[0] = 1.0 - w;
 			yWeight[0] *= yWeight[0];
 			yWeight[0] *= yWeight[0] * yWeight[0];
 			yWeight[0] *= (1.0 - w) / 5040.0;
 			w2 = w * w;
 			yWeight[1] = (120.0 / 7.0 + w * (-56.0 + w * (72.0 + w
 				* (-40.0 + w2 * (12.0 + w * (-6.0 + w)))))) / 720.0;
 			yWeight[2] = (397.0 / 7.0 - w * (245.0 / 3.0 + w * (-15.0 + w
 				* (-95.0 / 3.0 + w * (15.0 + w * (5.0 + w
 				* (-5.0 + w))))))) / 240.0;
 			yWeight[3] = (2416.0 / 35.0 + w2 * (-48.0 + w2 * (16.0 + w2
 				* (-4.0 + w)))) / 144.0;
 			yWeight[4] = (1191.0 / 35.0 - w * (-49.0 + w * (-9.0 + w
 				* (19.0 + w * (-3.0 + w) * (-3.0 + w2))))) / 144.0;
 			yWeight[5] = (40.0 / 7.0 + w * (56.0 / 3.0 + w * (24.0 + w
 				* (40.0 / 3.0 + w2 * (-4.0 + w * (-2.0 + w)))))) / 240.0;
 			yWeight[7] = w2;
 			yWeight[7] *= yWeight[7] * yWeight[7];
 			yWeight[7] *= w / 5040.0;
 			yWeight[6] = 1.0 - yWeight[0] - yWeight[1] - yWeight[2] - yWeight[3]
 				- yWeight[4] - yWeight[5] - yWeight[7];
             /* z */
 			w = coordinate[2] - (double)zIndex[3];
 			zWeight[0] = 1.0 - w;
 			zWeight[0] *= zWeight[0];
 			zWeight[0] *= zWeight[0] * zWeight[0];
 			zWeight[0] *= (1.0 - w) / 5040.0;
 			w2 = w * w;
 			zWeight[1] = (120.0 / 7.0 + w * (-56.0 + w * (72.0 + w
 				* (-40.0 + w2 * (12.0 + w * (-6.0 + w)))))) / 720.0;
 			zWeight[2] = (397.0 / 7.0 - w * (245.0 / 3.0 + w * (-15.0 + w
 				* (-95.0 / 3.0 + w * (15.0 + w * (5.0 + w
 				* (-5.0 + w))))))) / 240.0;
 			zWeight[3] = (2416.0 / 35.0 + w2 * (-48.0 + w2 * (16.0 + w2
 				* (-4.0 + w)))) / 144.0;
 			zWeight[4] = (1191.0 / 35.0 - w * (-49.0 + w * (-9.0 + w
 				* (19.0 + w * (-3.0 + w) * (-3.0 + w2))))) / 144.0;
 			zWeight[5] = (40.0 / 7.0 + w * (56.0 / 3.0 + w * (24.0 + w
 				* (40.0 / 3.0 + w2 * (-4.0 + w * (-2.0 + w)))))) / 240.0;
 			zWeight[7] = w2;
 			zWeight[7] *= zWeight[7] * zWeight[7];
 			zWeight[7] *= w / 5040.0;
 			zWeight[6] = 1.0 - zWeight[0] - zWeight[1] - zWeight[2] - zWeight[3]
 				- zWeight[4] - zWeight[5] - zWeight[7];
             /* t */
 			w = coordinate[3] - (double)tIndex[3];
 			tWeight[0] = 1.0 - w;
 			tWeight[0] *= tWeight[0];
 			tWeight[0] *= tWeight[0] * tWeight[0];
 			tWeight[0] *= (1.0 - w) / 5040.0;
 			w2 = w * w;
 			tWeight[1] = (120.0 / 7.0 + w * (-56.0 + w * (72.0 + w
 				* (-40.0 + w2 * (12.0 + w * (-6.0 + w)))))) / 720.0;
 			tWeight[2] = (397.0 / 7.0 - w * (245.0 / 3.0 + w * (-15.0 + w
 				* (-95.0 / 3.0 + w * (15.0 + w * (5.0 + w
 				* (-5.0 + w))))))) / 240.0;
 			tWeight[3] = (2416.0 / 35.0 + w2 * (-48.0 + w2 * (16.0 + w2
 				* (-4.0 + w)))) / 144.0;
 			tWeight[4] = (1191.0 / 35.0 - w * (-49.0 + w * (-9.0 + w
 				* (19.0 + w * (-3.0 + w) * (-3.0 + w2))))) / 144.0;
 			tWeight[5] = (40.0 / 7.0 + w * (56.0 / 3.0 + w * (24.0 + w
 				* (40.0 / 3.0 + w2 * (-4.0 + w * (-2.0 + w)))))) / 240.0;
 			tWeight[7] = w2;
 			tWeight[7] *= tWeight[7] * tWeight[7];
 			tWeight[7] *= w / 5040.0;
 			tWeight[6] = 1.0 - tWeight[0] - tWeight[1] - tWeight[2] - tWeight[3]
 				- tWeight[4] - tWeight[5] - tWeight[7];
 			break;
 		case 8L:
 			/* x */
 			w = coordinate[0] - (double)xIndex[4];
 			xWeight[0] = 1.0 / 2.0 - w;
 			xWeight[0] *= xWeight[0];
 			xWeight[0] *= xWeight[0];
 			xWeight[0] *= xWeight[0] / 40320.0;
 			w2 = w * w;
 			xWeight[1] = (39.0 / 16.0 - w * (6.0 + w * (-9.0 / 2.0 + w2)))
 				* (21.0 / 16.0 + w * (-15.0 / 4.0 + w * (9.0 / 2.0 + w
 				* (-3.0 + w)))) / 5040.0;
 			xWeight[2] = (82903.0 / 1792.0 + w * (-4177.0 / 32.0 + w
 				* (2275.0 / 16.0 + w * (-487.0 / 8.0 + w * (-85.0 / 8.0 + w
 				* (41.0 / 2.0 + w * (-5.0 + w * (-2.0 + w)))))))) / 1440.0;
 			xWeight[3] = (310661.0 / 1792.0 - w * (14219.0 / 64.0 + w
 				* (-199.0 / 8.0 + w * (-1327.0 / 16.0 + w * (245.0 / 8.0 + w
 				* (53.0 / 4.0 + w * (-8.0 + w * (-1.0 + w)))))))) / 720.0;
 			xWeight[4] = (2337507.0 / 8960.0 + w2 * (-2601.0 / 16.0 + w2
 				* (387.0 / 8.0 + w2 * (-9.0 + w2)))) / 576.0;
 			xWeight[5] = (310661.0 / 1792.0 - w * (-14219.0 / 64.0 + w
 				* (-199.0 / 8.0 + w * (1327.0 / 16.0 + w * (245.0 / 8.0 + w
 				* (-53.0 / 4.0 + w * (-8.0 + w * (1.0 + w)))))))) / 720.0;
 			xWeight[7] = (39.0 / 16.0 - w * (-6.0 + w * (-9.0 / 2.0 + w2)))
 				* (21.0 / 16.0 + w * (15.0 / 4.0 + w * (9.0 / 2.0 + w
 				* (3.0 + w)))) / 5040.0;
 			xWeight[8] = 1.0 / 2.0 + w;
 			xWeight[8] *= xWeight[8];
 			xWeight[8] *= xWeight[8];
 			xWeight[8] *= xWeight[8] / 40320.0;
 			xWeight[6] = 1.0 - xWeight[0] - xWeight[1] - xWeight[2] - xWeight[3]
 				- xWeight[4] - xWeight[5] - xWeight[7] - xWeight[8];
 			/* y */
 			w = coordinate[1] - (double)yIndex[4];
 			yWeight[0] = 1.0 / 2.0 - w;
 			yWeight[0] *= yWeight[0];
 			yWeight[0] *= yWeight[0];
 			yWeight[0] *= yWeight[0] / 40320.0;
 			w2 = w * w;
 			yWeight[1] = (39.0 / 16.0 - w * (6.0 + w * (-9.0 / 2.0 + w2)))
 				* (21.0 / 16.0 + w * (-15.0 / 4.0 + w * (9.0 / 2.0 + w
 				* (-3.0 + w)))) / 5040.0;
 			yWeight[2] = (82903.0 / 1792.0 + w * (-4177.0 / 32.0 + w
 				* (2275.0 / 16.0 + w * (-487.0 / 8.0 + w * (-85.0 / 8.0 + w
 				* (41.0 / 2.0 + w * (-5.0 + w * (-2.0 + w)))))))) / 1440.0;
 			yWeight[3] = (310661.0 / 1792.0 - w * (14219.0 / 64.0 + w
 				* (-199.0 / 8.0 + w * (-1327.0 / 16.0 + w * (245.0 / 8.0 + w
 				* (53.0 / 4.0 + w * (-8.0 + w * (-1.0 + w)))))))) / 720.0;
 			yWeight[4] = (2337507.0 / 8960.0 + w2 * (-2601.0 / 16.0 + w2
 				* (387.0 / 8.0 + w2 * (-9.0 + w2)))) / 576.0;
 			yWeight[5] = (310661.0 / 1792.0 - w * (-14219.0 / 64.0 + w
 				* (-199.0 / 8.0 + w * (1327.0 / 16.0 + w * (245.0 / 8.0 + w
 				* (-53.0 / 4.0 + w * (-8.0 + w * (1.0 + w)))))))) / 720.0;
 			yWeight[7] = (39.0 / 16.0 - w * (-6.0 + w * (-9.0 / 2.0 + w2)))
 				* (21.0 / 16.0 + w * (15.0 / 4.0 + w * (9.0 / 2.0 + w
 				* (3.0 + w)))) / 5040.0;
 			yWeight[8] = 1.0 / 2.0 + w;
 			yWeight[8] *= yWeight[8];
 			yWeight[8] *= yWeight[8];
 			yWeight[8] *= yWeight[8] / 40320.0;
 			yWeight[6] = 1.0 - yWeight[0] - yWeight[1] - yWeight[2] - yWeight[3]
 				- yWeight[4] - yWeight[5] - yWeight[7] - yWeight[8];
             /* z */
 			w = coordinate[2] - (double)zIndex[4];
 			zWeight[0] = 1.0 / 2.0 - w;
 			zWeight[0] *= zWeight[0];
 			zWeight[0] *= zWeight[0];
 			zWeight[0] *= zWeight[0] / 40320.0;
 			w2 = w * w;
 			zWeight[1] = (39.0 / 16.0 - w * (6.0 + w * (-9.0 / 2.0 + w2)))
 				* (21.0 / 16.0 + w * (-15.0 / 4.0 + w * (9.0 / 2.0 + w
 				* (-3.0 + w)))) / 5040.0;
 			zWeight[2] = (82903.0 / 1792.0 + w * (-4177.0 / 32.0 + w
 				* (2275.0 / 16.0 + w * (-487.0 / 8.0 + w * (-85.0 / 8.0 + w
 				* (41.0 / 2.0 + w * (-5.0 + w * (-2.0 + w)))))))) / 1440.0;
 			zWeight[3] = (310661.0 / 1792.0 - w * (14219.0 / 64.0 + w
 				* (-199.0 / 8.0 + w * (-1327.0 / 16.0 + w * (245.0 / 8.0 + w
 				* (53.0 / 4.0 + w * (-8.0 + w * (-1.0 + w)))))))) / 720.0;
 			zWeight[4] = (2337507.0 / 8960.0 + w2 * (-2601.0 / 16.0 + w2
 				* (387.0 / 8.0 + w2 * (-9.0 + w2)))) / 576.0;
 			zWeight[5] = (310661.0 / 1792.0 - w * (-14219.0 / 64.0 + w
 				* (-199.0 / 8.0 + w * (1327.0 / 16.0 + w * (245.0 / 8.0 + w
 				* (-53.0 / 4.0 + w * (-8.0 + w * (1.0 + w)))))))) / 720.0;
 			zWeight[7] = (39.0 / 16.0 - w * (-6.0 + w * (-9.0 / 2.0 + w2)))
 				* (21.0 / 16.0 + w * (15.0 / 4.0 + w * (9.0 / 2.0 + w
 				* (3.0 + w)))) / 5040.0;
 			zWeight[8] = 1.0 / 2.0 + w;
 			zWeight[8] *= zWeight[8];
 			zWeight[8] *= zWeight[8];
 			zWeight[8] *= zWeight[8] / 40320.0;
 			zWeight[6] = 1.0 - zWeight[0] - zWeight[1] - zWeight[2] - zWeight[3]
 				- zWeight[4] - zWeight[5] - zWeight[7] - zWeight[8];
             /* t */
 			w = coordinate[3] - (double)tIndex[4];
 			tWeight[0] = 1.0 / 2.0 - w;
 			tWeight[0] *= tWeight[0];
 			tWeight[0] *= tWeight[0];
 			tWeight[0] *= tWeight[0] / 40320.0;
 			w2 = w * w;
 			tWeight[1] = (39.0 / 16.0 - w * (6.0 + w * (-9.0 / 2.0 + w2)))
 				* (21.0 / 16.0 + w * (-15.0 / 4.0 + w * (9.0 / 2.0 + w
 				* (-3.0 + w)))) / 5040.0;
 			tWeight[2] = (82903.0 / 1792.0 + w * (-4177.0 / 32.0 + w
 				* (2275.0 / 16.0 + w * (-487.0 / 8.0 + w * (-85.0 / 8.0 + w
 				* (41.0 / 2.0 + w * (-5.0 + w * (-2.0 + w)))))))) / 1440.0;
 			tWeight[3] = (310661.0 / 1792.0 - w * (14219.0 / 64.0 + w
 				* (-199.0 / 8.0 + w * (-1327.0 / 16.0 + w * (245.0 / 8.0 + w
 				* (53.0 / 4.0 + w * (-8.0 + w * (-1.0 + w)))))))) / 720.0;
 			tWeight[4] = (2337507.0 / 8960.0 + w2 * (-2601.0 / 16.0 + w2
 				* (387.0 / 8.0 + w2 * (-9.0 + w2)))) / 576.0;
 			tWeight[5] = (310661.0 / 1792.0 - w * (-14219.0 / 64.0 + w
 				* (-199.0 / 8.0 + w * (1327.0 / 16.0 + w * (245.0 / 8.0 + w
 				* (-53.0 / 4.0 + w * (-8.0 + w * (1.0 + w)))))))) / 720.0;
 			tWeight[7] = (39.0 / 16.0 - w * (-6.0 + w * (-9.0 / 2.0 + w2)))
 				* (21.0 / 16.0 + w * (15.0 / 4.0 + w * (9.0 / 2.0 + w
 				* (3.0 + w)))) / 5040.0;
 			tWeight[8] = 1.0 / 2.0 + w;
 			tWeight[8] *= tWeight[8];
 			tWeight[8] *= tWeight[8];
 			tWeight[8] *= tWeight[8] / 40320.0;
 			tWeight[6] = 1.0 - tWeight[0] - tWeight[1] - tWeight[2] - tWeight[3]
 				- tWeight[4] - tWeight[5] - tWeight[7] - tWeight[8];
 			break;
 		case 9L:
 			/* x */
			w = coordinate[0] - (double)xIndex[4];						
			xWeight[0] = 1.0 - w;
			xWeight[0] *= xWeight[0];
			xWeight[0] *= xWeight[0];
			xWeight[0] *= xWeight[0] * (1.0 - w) / 362880.0;
			xWeight[1] = (502.0 / 9.0 + w * (-246.0 + w * (472.0 + w
				* (-504.0 + w * (308.0 + w * (-84.0 + w * (-56.0 / 3.0 + w
				* (24.0 + w * (-8.0 + w))))))))) / 40320.0;
			xWeight[2] = (3652.0 / 9.0 - w * (2023.0 / 2.0 + w * (-952.0 + w
				* (938.0 / 3.0 + w * (112.0 + w * (-119.0 + w * (56.0 / 3.0 + w
				* (14.0 + w * (-7.0 + w))))))))) / 10080.0;
			xWeight[3] = (44117.0 / 42.0 + w * (-2427.0 / 2.0 + w * (66.0 + w
				* (434.0 + w * (-129.0 + w * (-69.0 + w * (34.0 + w * (6.0 + w
				* (-6.0 + w))))))))) / 4320.0;
			w2 = w * w;
			xWeight[4] = (78095.0 / 63.0 - w2 * (700.0 + w2 * (-190.0 + w2
				* (100.0 / 3.0 + w2 * (-5.0 + w))))) / 2880.0;
			xWeight[5] = (44117.0 / 63.0 + w * (809.0 + w * (44.0 + w
				* (-868.0 / 3.0 + w * (-86.0 + w * (46.0 + w * (68.0 / 3.0 + w
				* (-4.0 + w * (-4.0 + w))))))))) / 2880.0;
			xWeight[6] = (3652.0 / 21.0 - w * (-867.0 / 2.0 + w * (-408.0 + w
				* (-134.0 + w * (48.0 + w * (51.0 + w * (-4.0 + w) * (-1.0 + w)
				* (2.0 + w))))))) / 4320.0;
			xWeight[7] = (251.0 / 18.0 + w * (123.0 / 2.0 + w * (118.0 + w
				* (126.0 + w * (77.0 + w * (21.0 + w * (-14.0 / 3.0 + w
				* (-6.0 + w * (-2.0 + w))))))))) / 10080.0;
			xWeight[9] = w2 * w2;
			xWeight[9] *= xWeight[9] * w / 362880.0;
			xWeight[8] = 1.0 - xWeight[0] - xWeight[1] - xWeight[2] - xWeight[3]
				- xWeight[4] - xWeight[5] - xWeight[6] - xWeight[7] - xWeight[9];
			/* y */
			w = coordinate[1] - (double)yIndex[4];
			yWeight[0] = 1.0 - w;
			yWeight[0] *= yWeight[0];
			yWeight[0] *= yWeight[0];
			yWeight[0] *= yWeight[0] * (1.0 - w) / 362880.0;
			yWeight[1] = (502.0 / 9.0 + w * (-246.0 + w * (472.0 + w
				* (-504.0 + w * (308.0 + w * (-84.0 + w * (-56.0 / 3.0 + w
				* (24.0 + w * (-8.0 + w))))))))) / 40320.0;
			yWeight[2] = (3652.0 / 9.0 - w * (2023.0 / 2.0 + w * (-952.0 + w
				* (938.0 / 3.0 + w * (112.0 + w * (-119.0 + w * (56.0 / 3.0 + w
				* (14.0 + w * (-7.0 + w))))))))) / 10080.0;
			yWeight[3] = (44117.0 / 42.0 + w * (-2427.0 / 2.0 + w * (66.0 + w
				* (434.0 + w * (-129.0 + w * (-69.0 + w * (34.0 + w * (6.0 + w
				* (-6.0 + w))))))))) / 4320.0;
			w2 = w * w;
			yWeight[4] = (78095.0 / 63.0 - w2 * (700.0 + w2 * (-190.0 + w2
				* (100.0 / 3.0 + w2 * (-5.0 + w))))) / 2880.0;
			yWeight[5] = (44117.0 / 63.0 + w * (809.0 + w * (44.0 + w
				* (-868.0 / 3.0 + w * (-86.0 + w * (46.0 + w * (68.0 / 3.0 + w
				* (-4.0 + w * (-4.0 + w))))))))) / 2880.0;
			yWeight[6] = (3652.0 / 21.0 - w * (-867.0 / 2.0 + w * (-408.0 + w
				* (-134.0 + w * (48.0 + w * (51.0 + w * (-4.0 + w) * (-1.0 + w)
				* (2.0 + w))))))) / 4320.0;
			yWeight[7] = (251.0 / 18.0 + w * (123.0 / 2.0 + w * (118.0 + w
				* (126.0 + w * (77.0 + w * (21.0 + w * (-14.0 / 3.0 + w
				* (-6.0 + w * (-2.0 + w))))))))) / 10080.0;
			yWeight[9] = w2 * w2;
			yWeight[9] *= yWeight[9] * w / 362880.0;
			yWeight[8] = 1.0 - yWeight[0] - yWeight[1] - yWeight[2] - yWeight[3]
				- yWeight[4] - yWeight[5] - yWeight[6] - yWeight[7] - yWeight[9];
             /* z */
 			w = coordinate[2] - (double)zIndex[4];
			zWeight[0] = 1.0 - w;
			zWeight[0] *= zWeight[0];
			zWeight[0] *= zWeight[0];
			zWeight[0] *= zWeight[0] * (1.0 - w) / 362880.0;
			zWeight[1] = (502.0 / 9.0 + w * (-246.0 + w * (472.0 + w
				* (-504.0 + w * (308.0 + w * (-84.0 + w * (-56.0 / 3.0 + w
				* (24.0 + w * (-8.0 + w))))))))) / 40320.0;
			zWeight[2] = (3652.0 / 9.0 - w * (2023.0 / 2.0 + w * (-952.0 + w
				* (938.0 / 3.0 + w * (112.0 + w * (-119.0 + w * (56.0 / 3.0 + w
				* (14.0 + w * (-7.0 + w))))))))) / 10080.0;
			zWeight[3] = (44117.0 / 42.0 + w * (-2427.0 / 2.0 + w * (66.0 + w
				* (434.0 + w * (-129.0 + w * (-69.0 + w * (34.0 + w * (6.0 + w
				* (-6.0 + w))))))))) / 4320.0;
			w2 = w * w;
			zWeight[4] = (78095.0 / 63.0 - w2 * (700.0 + w2 * (-190.0 + w2
				* (100.0 / 3.0 + w2 * (-5.0 + w))))) / 2880.0;
			zWeight[5] = (44117.0 / 63.0 + w * (809.0 + w * (44.0 + w
				* (-868.0 / 3.0 + w * (-86.0 + w * (46.0 + w * (68.0 / 3.0 + w
				* (-4.0 + w * (-4.0 + w))))))))) / 2880.0;
			zWeight[6] = (3652.0 / 21.0 - w * (-867.0 / 2.0 + w * (-408.0 + w
				* (-134.0 + w * (48.0 + w * (51.0 + w * (-4.0 + w) * (-1.0 + w)
				* (2.0 + w))))))) / 4320.0;
			zWeight[7] = (251.0 / 18.0 + w * (123.0 / 2.0 + w * (118.0 + w
				* (126.0 + w * (77.0 + w * (21.0 + w * (-14.0 / 3.0 + w
				* (-6.0 + w * (-2.0 + w))))))))) / 10080.0;
			zWeight[9] = w2 * w2;
			zWeight[9] *= zWeight[9] * w / 362880.0;
			zWeight[8] = 1.0 - zWeight[0] - zWeight[1] - zWeight[2] - zWeight[3]
				- zWeight[4] - zWeight[5] - zWeight[6] - zWeight[7] - zWeight[9];
             /* t */
 			w = coordinate[3] - (double)tIndex[4];
			tWeight[0] = 1.0 - w;
			tWeight[0] *= tWeight[0];
			tWeight[0] *= tWeight[0];
			tWeight[0] *= tWeight[0] * (1.0 - w) / 362880.0;
			tWeight[1] = (502.0 / 9.0 + w * (-246.0 + w * (472.0 + w
				* (-504.0 + w * (308.0 + w * (-84.0 + w * (-56.0 / 3.0 + w
				* (24.0 + w * (-8.0 + w))))))))) / 40320.0;
			tWeight[2] = (3652.0 / 9.0 - w * (2023.0 / 2.0 + w * (-952.0 + w
				* (938.0 / 3.0 + w * (112.0 + w * (-119.0 + w * (56.0 / 3.0 + w
				* (14.0 + w * (-7.0 + w))))))))) / 10080.0;
			tWeight[3] = (44117.0 / 42.0 + w * (-2427.0 / 2.0 + w * (66.0 + w
				* (434.0 + w * (-129.0 + w * (-69.0 + w * (34.0 + w * (6.0 + w
				* (-6.0 + w))))))))) / 4320.0;
			w2 = w * w;
			tWeight[4] = (78095.0 / 63.0 - w2 * (700.0 + w2 * (-190.0 + w2
				* (100.0 / 3.0 + w2 * (-5.0 + w))))) / 2880.0;
			tWeight[5] = (44117.0 / 63.0 + w * (809.0 + w * (44.0 + w
				* (-868.0 / 3.0 + w * (-86.0 + w * (46.0 + w * (68.0 / 3.0 + w
				* (-4.0 + w * (-4.0 + w))))))))) / 2880.0;
			tWeight[6] = (3652.0 / 21.0 - w * (-867.0 / 2.0 + w * (-408.0 + w
				* (-134.0 + w * (48.0 + w * (51.0 + w * (-4.0 + w) * (-1.0 + w)
				* (2.0 + w))))))) / 4320.0;
			tWeight[7] = (251.0 / 18.0 + w * (123.0 / 2.0 + w * (118.0 + w
				* (126.0 + w * (77.0 + w * (21.0 + w * (-14.0 / 3.0 + w
				* (-6.0 + w * (-2.0 + w))))))))) / 10080.0;
			tWeight[9] = w2 * w2;
			tWeight[9] *= tWeight[9] * w / 362880.0;
			tWeight[8] = 1.0 - tWeight[0] - tWeight[1] - tWeight[2] - tWeight[3]
				- tWeight[4] - tWeight[5] - tWeight[6] - tWeight[7] - tWeight[9];
 			break;

 		case -3L: // oMom3
 			/* x */
                xWeight[0] = oMom3(coordinate[0] - (double)xIndex[0]); 
                xWeight[1] = oMom3(coordinate[0] - (double)xIndex[1]);
                xWeight[2] = oMom3(coordinate[0] - (double)xIndex[2]);
                xWeight[3] = oMom3(coordinate[0] - (double)xIndex[3]);
				/* y */
                yWeight[0] = oMom3(coordinate[1] - (double)yIndex[0]); 
                yWeight[1] = oMom3(coordinate[1] - (double)yIndex[1]);
                yWeight[2] = oMom3(coordinate[1] - (double)yIndex[2]);
                yWeight[3] = oMom3(coordinate[1] - (double)yIndex[3]);
				/* z */
                zWeight[0] = oMom3(coordinate[2] - (double)zIndex[0]); 
                zWeight[1] = oMom3(coordinate[2] - (double)zIndex[1]);
                zWeight[2] = oMom3(coordinate[2] - (double)zIndex[2]);
                zWeight[3] = oMom3(coordinate[2] - (double)zIndex[3]);
				/* t */
                tWeight[0] = oMom3(coordinate[3] - (double)tIndex[0]); 
                tWeight[1] = oMom3(coordinate[3] - (double)tIndex[1]);
                tWeight[2] = oMom3(coordinate[3] - (double)tIndex[2]);
                tWeight[3] = oMom3(coordinate[3] - (double)tIndex[3]);
 			break;

 		default:
 			printf("Invalid spline degree\n");
 			break;
 	}

    //printf("Size4 = %ld %ld %ld %ld\n",Size[0],Size[1],Size[2],Size[3]);
     //printf("indices[PRE] = [%d,%d]\n",xIndex[0],xIndex[1]);
     //printf("size_X  size_X2= [%d,%d]\n",size_X,size_X2);
 	/* apply the mirror boundary conditions */
    
 	for (k = 0L; k <= abs(SplineDegree); k++) {
 		xIndex[k] = (size_X == 1L) ? (0L) : ((xIndex[k] < 0L) ?
 			(-xIndex[k] - size_X2 * ((-xIndex[k]) / size_X2))
 			: (xIndex[k] - size_X2 * (xIndex[k] / size_X2)));
 		if (size_X <= xIndex[k]) {
 			xIndex[k] = size_X2 - xIndex[k];
 		}
 		yIndex[k] = (size_Y == 1L) ? (0L) : ((yIndex[k] < 0L) ?
 			(-yIndex[k] - size_Y2 * ((-yIndex[k]) / size_Y2))
 			: (yIndex[k] - size_Y2 * (yIndex[k] / size_Y2)));
 		if (size_Y <= yIndex[k]) {
 			yIndex[k] = size_Y2 - yIndex[k];
 		}
         zIndex[k] = (size_Z == 1L) ? (0L) : ((zIndex[k] < 0L) ?
 			(-zIndex[k] - size_Z2 * ((-zIndex[k]) / size_Z2))
 			: (zIndex[k] - size_Z2 * (zIndex[k] / size_Z2)));
 		if (size_Z <= zIndex[k]) {
 			zIndex[k] = size_Z2 - zIndex[k];
 		}
 		tIndex[k] = (size_T == 1L) ? (0L) : ((tIndex[k] < 0L) ?
 			(-tIndex[k] - size_T2 * ((-tIndex[k]) / size_T2))
 			: (tIndex[k] - size_T2 * (tIndex[k] / size_T2)));
 		if (size_T <= tIndex[k]) {
 			tIndex[k] = size_T2 - tIndex[k];
 		}

 	}
     //printf("Size5 = %ld %ld %ld %ld\n",Size[0],Size[1],Size[2],Size[3]);
    // printf("indices[POST] = [%d,%d]\n",xIndex[0],xIndex[1]);

    
//     printf("xIndex[] = [%ld,%ld,%ld,%ld]\n",xIndex[0],xIndex[1],xIndex[2],xIndex[3]);
//     printf("xWeight[] = [%f,%f,%f,%f]\n",xWeight[0],xWeight[1],xWeight[2],xWeight[3]);
//     printf("yIndex[] = [%ld,%ld,%ld,%ld]\n",yIndex[0],yIndex[1],yIndex[2],yIndex[3]);
//     printf("yWeight[] = [%f,%f,%f,%f]\n",yWeight[0],yWeight[1],yWeight[2],xWeight[3]);
//     printf("zIndex[] = [%ld,%ld,%ld,%ld]\n",zIndex[0],zIndex[1],zIndex[2],zIndex[3]);
//     printf("zWeight[] = [%f,%f,%f,%f]\n",zWeight[0],zWeight[1],zWeight[2],zWeight[3]);
//     printf("tIndex[] = [%ld,%ld,%ld,%ld]\n",tIndex[0],tIndex[1],tIndex[2],tIndex[3]);
//     printf("tWeight[] = [%f,%f,%f,%f]\n",tWeight[0],tWeight[1],tWeight[2],tWeight[3]);
    
 	/* perform interpolation */
 	interpolated = 0.0;


      for(n = 0L; n <= abs(SplineDegree); n++){
          w = 0.0;
          w2 = 0.0;
          w3 = 0.0;
          for(k = 0L; k<= abs(SplineDegree); k++){
             for(j = 0L; j<= abs(SplineDegree); j++){
                 for (i = 0L; i <= abs(SplineDegree); i++) {
                      w += xWeight[i] * Bcoeff[xIndex[i] + (yIndex[j]*size_X) + (zIndex[k]*size_X*size_Y) + (tIndex[n]*size_X*size_Y*size_Z)];
                 }
                 w2 += yWeight[j] * w;
                 w = 0.0;
             }
             w3 += zWeight[k] * w2;
             w2 = 0.0;
          }
          interpolated += tWeight[n] * w3;
          w3 = 0.0;
     }
//     
//     printf("coordenadaX[] = [%f]\n",coordinate[0]);
//             printf("\tindicesX[] = [%ld,%ld]\n",xIndex[0],xIndex[1]);
//             printf("\tpesosX[] = [%f,%f]\n",xWeight[0],xWeight[1]);
//     printf("coordenaday[] = [%f]\n",coordinate[1]);
//             printf("\tindicesY[] = [%ld,%ld]\n",yIndex[0],yIndex[1]);
//             printf("\tpesosY[] = [%f,%f]\n",yWeight[0],yWeight[1]);
//             printf("\tvaloresY[0,0] = [%f]\n",Bcoeff[xIndex[0] + (yIndex[0]*size_X)]);
//             printf("\tvaloresY[0,1] = [%f]\n",Bcoeff[xIndex[0] + (yIndex[1]*size_X)]);
//             printf("\tvaloresY[1,0] = [%f]\n",Bcoeff[xIndex[1] + (yIndex[0]*size_X)]);
//             printf("\tvaloresY[1,1] = [%f]\n",Bcoeff[xIndex[1] + (yIndex[1]*size_X)]);
//     printf("Interpolated = %f\n",interpolated);
//     
//     
    
      return interpolated;

}

extern double oMom3(
					double x				/*input for the oMom3 kernel */
				){

    //printf("Valor de x antes = %f\n",(x));
	x = (x<0.0)?-x:x;

	if(x < 1.0){

		return (0.5*x*x*x - x*x + x/14.0 +13.0/21.0);

	}else if (x < 2.0){

		return (- x*x*x/6 + x*x - x*85.0/42.0 + 29.0/21.0);

		}else{
			return 0;
		}

}

extern double lin(
					double x				/*input for the lineal kernel */
				){
    
    double x_abs = (x>0.0)?x:(-x);
//     printf("Dato entrada %f\n",x);
//     printf("Dato entrada abs %f\n",x_abs);
    
    if(x_abs>1.0){
        return 0;
    }else{
        return (1.0-x_abs);
    }
    
}















