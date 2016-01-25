#include "grid.h"

void nn_interp_to_mesh(e *E, int n, double weight, point *pin, int nx, int ny, point *interp){


	double	xmin = NaN;
	double	xmax = NaN;
 	double	ymin = NaN;
	double 	ymax = NaN;
	point	*pout = NULL;
	delaunay *d = NULL;			// for the triangle functions
	void *interpolator = NULL;	// the interpolation method
	preader *pr = NULL;			// point reader struct
	int	this_point = 0;

	nn_rule = SIBSON;

	// get the range of the input data...
	if( (E->have_min_lat == TRUE) && (E->have_max_lat == TRUE)
		&&
		(E->have_min_lon == TRUE) && (E->have_max_lon == TRUE)){
			xmin = E->min_lon;
			xmax = E->max_lon;
			ymin = E->min_lat;
			ymax = E->max_lat;
	}
	else
		points_getrange(n, pin, 1.0, &xmin, &xmax, &ymin, &ymax);


	// creates the output regular grid x, y coordinates
    pr = preader_create1(xmin, xmax, ymin, ymax, nx, ny);
	// generate the mesh

	d = delaunay_build(n, pin, 0, NULL, 0, NULL);	// construct the triangulation

	interpolator = nnpi_create(d);
    //nnpi_setwmin(interpolator, -DBL_MAX);
	nnpi_setwmin(interpolator, weight);

	this_point = 0;
 	while ((pout = preader_getpoint(pr)) != NULL) {

		nnpi_interpolate_point(interpolator, pout);

		interp[this_point].x = pout->x;
		interp[this_point].y = pout->y;
		interp[this_point].z = pout->z;

		this_point++;

		//points_write(1, pout);

    }


	free(pout);
	free(d);			// for the triangle functions
	free(interpolator);	// the interpolation method
	free(pr);

}

void get_mesh_dimensions(e *E, int n, point *pin){

	double	xmin = NaN;
	double	xmax = NaN;
 	double	ymin = NaN;
	double 	ymax = NaN;


	// get the range of the input data...
	if( (E->have_min_lat == TRUE) && (E->have_max_lat == TRUE)
		&&
		(E->have_min_lon == TRUE) && (E->have_max_lon == TRUE)){
			xmin = E->min_lon;
			xmax = E->max_lon;
			ymin = E->min_lat;
			ymax = E->max_lat;
	}
	else
		points_getrange(n, pin, 1.0, &xmin, &xmax, &ymin, &ymax);


	E->nn_nx = ceil(1.0+(xmax - xmin)/E->nn_dx);
	E->nn_ny = ceil(1.0+(ymax - ymin)/E->nn_dy);

}
