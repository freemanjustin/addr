/******************************************************************************
 *
 * File:           minell.h
 *
 * Created:        24/02/2003
 *
 * Author:         Pavel Sakov
 *                 CSIRO Marine Research
 *
 * Purpose:        A header for the minimal ellipse stuff
 *
 * Revisions:      None
 *
 *****************************************************************************/

#if !defined(_MINELL_H)
#define _MINELL_H

#if !defined(_POINT_STRUCT)
#define _POINT_STRUCT
typedef struct point {
    double x;
    double y;
    double z;
	double ave;
	struct point *next;
} point;
#endif

#if !defined(_MINELL_STRUCT)
#define _MINELL_STRUCT
struct minell;
typedef struct minell minell;
#endif

/* Note that minell_build() shuffles the input point array */
minell* minell_build(int n, point p[]);
void minell_destroy(minell* me);
void minell_scalepoints(minell* me, int n, point p[]);
void minell_rescalepoints(minell* me, int n, point p[]);

#endif
