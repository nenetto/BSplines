
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
 * Laboratorio de Imagen M�dica
 * WEB del LIM
 * Unidad de Medicina y Cirug�a Experimental.
 * Fundaci�n para la Investigaci�n Biom�dica del
 * Hospital Gregorio Mara��n
 * C/ Doctor Esquerdo, 46
 * 28007  Madrid
 * Tel�fono: +34 91 4265017
 * Fax: +34 91 4265108
 *
****************************************************************************/




/*****************************************************************************
 *	System includes
 ****************************************************************************/
#include	<float.h>
#include	<math.h>
#include	<stddef.h>
#include	<stdio.h>
#include	<stdlib.h>
#include    "mex.h"
#include	"bspline_coeff.h"

/*****************************************************************************
 *	Define
 ****************************************************************************/

#define TOL 2.22044604925031308085e-16
