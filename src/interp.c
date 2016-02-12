// gridr
//
// freeman.justin@gmail.com


#include "grid.h"

#define TOLERANCE	1e-6
#define ABS(x)		((x) < 0 ? -(x) : (x))
#define MAX(a,b)    ( ((a) > (b)) ? (a) : (b) )
#define MIN(a,b)    ( ((a) < (b)) ? (a) : (b) )

double evaluate_linear_quad_shape_function( double *xi, double *eta, double *pos_to_interp_at, int node ){

	return ( 0.25 * (1.0 + pos_to_interp_at[0] * xi[node]) * (1.0 + pos_to_interp_at[1] * eta[node]) );

}

double relative_difference(double a, double b)
{
	double c = ABS(a);
	double d = ABS(b);

	d = MAX(c, d);

	return d == 0.0 ? 0.0 : ABS(a - b) / d;
}


void init_xi_eta(e *E){

	// anti-clockwise element numbering
	// master element nodal positions
    //
    //  E->xi is the x coordinate
    //  E->eta is the y coordinate
    //
    //            ^
    // (-1,1)     |      (1,1)
    //    o----------------o
    //    |       |        |
    //    |     (0,0)------|- >
    //    |                |
    //    o----------------o
    // (-1,-1)           (1,-1)
    //

	E->xi[0] =	-1.0;
	E->xi[1] =	 1.0;
	E->xi[2] =	 1.0;
	E->xi[3] =	-1.0;

	E->eta[0] = -1.0;
	E->eta[1] = -1.0;
	E->eta[2] =  1.0;
	E->eta[3] =  1.0;

}

void calculate_interpolation_weights(element *el, double *xi, double *eta, double *pos){

	int i;
	double	pos_local[2];

	// transform from global interpolation position to local element coordinates
	pos_local[0] = (pos[0] - 0.5*(el->node_coord[1][0] + el->node_coord[0][0]))/( 0.5 * (el->node_coord[1][0]-el->node_coord[0][0]) );
	pos_local[1] = (pos[1] - 0.5*(el->node_coord[3][1] + el->node_coord[0][1]))/( 0.5 * (el->node_coord[3][1]-el->node_coord[0][1]) );

	for(i=0; i < 4; i++)
		el->interp_weights[i] = evaluate_linear_quad_shape_function( xi, eta, pos_local, i );

}


void interpolate_point(element *el, double *interp_value){

	int i;

	*interp_value = 0.0;
	for(i=0; i < 4; i++) {
		*interp_value += el->node_value[i] * el->interp_weights[i];
    }
}



void store_mesh(e *E, int nx, int ny){

	int i,el;

	el=0;

	for(i=0;i<nx;i++){
		E->msh.x[i] = E->ele[i].node_coord[0][0];
	}
	E->msh.x[nx] = E->ele[nx-1].node_coord[1][0];

	for(i=0;i<ny;i++){
		E->msh.y[i] = E->ele[i*nx].node_coord[0][1];
	}
	E->msh.y[ny] = E->ele[(ny-1)*nx].node_coord[3][1];
}

int	get_owner_element(e *E, double *pos){

	int		i;
	int		max_its;
	int		element;
	int		element_new;
	int		x_element;
	int		y_element;
	double	point;
	double	interval_size;
	int		check = 1;


	max_its = 20;

	// search over elements to determine which element contains this position
	// apply a binary search method on the elements

	// search in x
    i=0;
    element = (double)E->nx/2.0;
    element_new = 0;
	interval_size = E->nx/2.0;

	do{
		if( (relative_difference(pos[0],E->msh.x[element+1])<=TOLERANCE) || (relative_difference(pos[0],E->msh.x[element])<=TOLERANCE)
		   ||
		   (pos[0] <= E->msh.x[element+1]) && (pos[0] >= E->msh.x[element]) ){
			// found it
			element = element;
			break;
		}
		else if( E->msh.x[element] >= pos[0] ){ // means that the pos point is on the left of our current index
			element_new = ( (double)element - ceil(interval_size/2.0) );
			if(element_new < 0) element_new=0;
		}
		else{	// means that the pos point is on the right of our current index

			element_new = ceil( (double)element + (interval_size/2.0) );

			if(element_new >= E->nx){
				element_new=E->nx-1;
			}
		}
		interval_size = fabs(element_new - element);
		element = element_new;
		i++;
	}while(i<max_its);

	if(i>=max_its){
		printf("x get owner failed!\n");
		printf("pos = %f, %f\n", pos[0], pos[1]);
		printf("its = %d\n", i);

		printf(" (%.5f, %.5f) --------- (%.5f, %.5f)\n", E->ele[element].node_coord[3][0], E->ele[element].node_coord[3][1],
			   E->ele[element].node_coord[2][0], E->ele[element].node_coord[2][1]);
		printf("   |                              |\n");
		printf("   |                              |\n");
		printf("   |                              |\n");
		printf("   |                              |\n");
		printf(" (%.5f, %.5f) --------- (%.5f, %.5f)\n\n", E->ele[element].node_coord[0][0], E->ele[element].node_coord[0][1],
			   E->ele[element].node_coord[1][0], E->ele[element].node_coord[1][1]);

		exit(1);
	}

	// x_element stores the column that this point will be in
	// lets use this information to constrain our search in y

	x_element = element;
	element = ceil((double)E->ny/2.0);
	element_new = 0;//(double)E->ny/2.0;
	interval_size = E->ny/2.0;

	// now search in y
	i=0;
	do{
		if( (relative_difference(pos[1],E->msh.y[element+1])<=TOLERANCE) || (relative_difference(pos[1],E->msh.y[element])<=TOLERANCE)
		   ||
		   (pos[1] <= E->msh.y[element+1]) && (pos[1] >= E->msh.y[element]) ){
			// found it
			break;
		}

		if( (E->msh.y[element] >= pos[1]) ){
			element_new = ( (double)element - ceil(interval_size/2.0) ) ;
			if(element_new < 0) element_new=0;
		}
		else{
			element_new = ( (double)element + ceil(interval_size/2.0) );

			if(element_new >= E->ny) {
				element_new=E->ny-1;
			}
        }

		interval_size = fabs(element_new - element);
		element = element_new;


		i++;
	}while(i<max_its);

	if(i>=max_its){
		printf("y get owner failed!\n");
		printf("pos = %f, %f\n", pos[0], pos[1]);
		printf("its = %d\n", i);

		printf(" (%.5f, %.5f) --------- (%.5f, %.5f)\n", E->ele[element].node_coord[3][0], E->ele[element].node_coord[3][1],
			   E->ele[element].node_coord[2][0], E->ele[element].node_coord[2][1]);
		printf("   |                              |\n");
		printf("   |                              |\n");
		printf("   |                              |\n");
		printf("   |                              |\n");
		printf(" (%.5f, %.5f) --------- (%.5f, %.5f)\n\n", E->ele[element].node_coord[0][0], E->ele[element].node_coord[0][1],
			   E->ele[element].node_coord[1][0], E->ele[element].node_coord[1][1]);

		exit(1);
	}

	// the owner element for this point is element * x_element
	element = (element*(E->nx)) + x_element;

    /*
        printf("get owner finished\n");
		printf("pos = %f, %f\n", pos[0], pos[1]);
		printf("its = %d\n", i);

		printf(" (%.5f, %.5f) --------- (%.5f, %.5f)\n", E->ele[element].node_coord[3][0], E->ele[element].node_coord[3][1],
			   E->ele[element].node_coord[2][0], E->ele[element].node_coord[2][1]);
		printf("   |                              |\n");
		printf("   |                              |\n");
		printf("   |                              |\n");
		printf("   |                              |\n");
		printf(" (%.5f, %.5f) --------- (%.5f, %.5f)\n\n", E->ele[element].node_coord[0][0], E->ele[element].node_coord[0][1],
			   E->ele[element].node_coord[1][0], E->ele[element].node_coord[1][1]);
    */

	return element;
}

void print_elements(e *E){

	int i;

	for(i=0;i<E->nElements;i++){
		printf("Element %d:\n", i);
		printf("\tnode pos: 0 = (%f,%f) 1 = (%f,%f) 2 = (%f,%f) 3 = (%f,%f)\n", E->ele[i].node_coord[0][0], E->ele[i].node_coord[0][1],E->ele[i].node_coord[1][0], E->ele[i].node_coord[1][1],
               E->ele[i].node_coord[2][0], E->ele[i].node_coord[2][1],E->ele[i].node_coord[3][0], E->ele[i].node_coord[3][1]);
		printf("\tnode values: 0 = %f, 1 = %f, 2 = %f, 3 = %f\n", E->ele[i].node_value[0], E->ele[i].node_value[1], E->ele[i].node_value[2], E->ele[i].node_value[3]);

	}

}



void interp_bathy_on_grid(e *E){

    int i,j;
    int el;
    double  pos[2];
    //double  interp_value;


    // setup the interpolation source grid data structures
	E->nx = E->b.nlon-1; //E->nc.x-1;
	E->ny = E->b.nlat-1; //E->nc.y-1;
	E->nElements = E->nx*E->ny ;
	E->ele = malloc(E->nElements*sizeof(element));
	E->nodesPerEl = 4;

	E->msh.x = malloc((E->nx+1) * sizeof(double));
	E->msh.y = malloc((E->ny+1) * sizeof(double));

	init_xi_eta(E);

	el = 0;
	for(i=0;i<E->ny;i++){
		for(j=0;j<E->nx;j++){

			// set positions for this element
			// element node numbering is anti-clockwise
			E->ele[el].node_coord[0][0] = E->b.lon[j];
			E->ele[el].node_coord[0][1] = E->b.lat[i]; // 0

			E->ele[el].node_coord[1][0] = E->b.lon[j+1];
			E->ele[el].node_coord[1][1] = E->b.lat[i]; // 1

			E->ele[el].node_coord[2][0] = E->b.lon[j+1];
			E->ele[el].node_coord[2][1] = E->b.lat[i+1]; // 2

			E->ele[el].node_coord[3][0] = E->b.lon[j];
			E->ele[el].node_coord[3][1] = E->b.lat[i+1]; // 3

			// now set the nodal value for this element
			E->ele[el].node_value[0] = -E->b.field[i][j];
			E->ele[el].node_value[1] = -E->b.field[i][j+1];
			E->ele[el].node_value[2] = -E->b.field[i+1][j+1];
			E->ele[el].node_value[3] = -E->b.field[i+1][j];

			el++;
		}
	}

    // the store_mesh call is required for the function get_owner_element
	store_mesh(E, E->nx, E->ny);
    //print_elements(E);

    // for each rho grid point
    for(i=0;i<E->g.nY;i++){
        for(j=0;j<E->g.nX;j++){

            // get pos of grid point
            pos[0] = E->lon_rho[i][j];
            pos[1] = E->lat_rho[i][j];

            // find out which element this lies within
			el = get_owner_element(E, pos);

            calculate_interpolation_weights(&E->ele[el], E->xi, E->eta, pos);
			interpolate_point(&E->ele[el], &E->h[i][j]);

        }
    }
}



void interp_tide_to_roms(e *E, int t){

    int i,j;
    int el;
    double  pos[2];
    //double  interp_value;


    // setup the interpolation source grid data structures
	// source grid is the tide data
	E->nx = E->nLonTide-1; //E->nc.x-1;
	E->ny = E->nLatTide-1; //E->nc.y-1;
	E->nElements = E->nx*E->ny ;
	E->ele = malloc(E->nElements*sizeof(element));
	E->nodesPerEl = 4;

	E->msh.x = malloc((E->nx+1) * sizeof(double));
	E->msh.y = malloc((E->ny+1) * sizeof(double));

	init_xi_eta(E);

	el = 0;
	for(i=0;i<E->ny;i++){
		for(j=0;j<E->nx;j++){

			// set positions for this element
			// element node numbering is anti-clockwise
			E->ele[el].node_coord[0][0] = E->tideLon[j];
			E->ele[el].node_coord[0][1] = E->tideLat[i]; // 0

			E->ele[el].node_coord[1][0] = E->tideLon[j+1];
			E->ele[el].node_coord[1][1] = E->tideLat[i]; // 1

			E->ele[el].node_coord[2][0] = E->tideLon[j+1];
			E->ele[el].node_coord[2][1] = E->tideLat[i+1]; // 2

			E->ele[el].node_coord[3][0] = E->tideLon[j];
			E->ele[el].node_coord[3][1] = E->tideLat[i+1]; // 3

			// now set the nodal value for this element
			E->ele[el].node_value[0] = E->tide_data[t][0][i][j];
			E->ele[el].node_value[1] = E->tide_data[t][0][i][j+1];
			E->ele[el].node_value[2] = E->tide_data[t][0][i+1][j+1];
			E->ele[el].node_value[3] = E->tide_data[t][0][i+1][j];

			el++;
		}
	}

    // the store_mesh call is required for the function get_owner_element
	store_mesh(E, E->nx, E->ny);
    //print_elements(E);

    // for each rho grid point
    for(i=0;i<E->nLonRho;i++){
        for(j=0;j<E->nLatRho;j++){
            // get pos of grid point
            pos[0] = E->lon_rho[i][j];
            pos[1] = E->lat_rho[i][j];

			/*
			if(pos[0] > 156.0){
				//printf("pre: i = %d, j = %d, pos = %f, %f\n",i,j,pos[0],pos[1]);
				pos[0] = 156.0;
				//printf("  modded: i = %d, j = %d, pos = %f, %f\n",i,j,pos[0],pos[1]);
			}

			if(pos[1] > -9.0){
				//printf("pre: i = %d, j = %d, pos = %f, %f\n",i,j,pos[0],pos[1]);
				pos[1] = -9.0;
				//printf("  modded: i = %d, j = %d, pos = %f, %f\n",i,j,pos[0],pos[1]);
			}
			*/

            // find out which element this lies within
			el = get_owner_element(E, pos);

            calculate_interpolation_weights(&E->ele[el], E->xi, E->eta, pos);
			interpolate_point(&E->ele[el], &E->tide_on_roms[t][i][j]);

        }
    }
}


void interp_hs_to_roms(e *E, int t){

    int i,j;
    int el;
    double  pos[2];
    //double  interp_value;


    // setup the interpolation source grid data structures
	// source grid is the tide data
	E->nx = E->nLonWaves-1; //E->nc.x-1;
	E->ny = E->nLatWaves-1; //E->nc.y-1;
	E->nElements = E->nx*E->ny ;
	E->ele = malloc(E->nElements*sizeof(element));
	E->nodesPerEl = 4;

	E->msh.x = malloc((E->nx+1) * sizeof(double));
	E->msh.y = malloc((E->ny+1) * sizeof(double));

	init_xi_eta(E);

	el = 0;
	for(i=0;i<E->ny;i++){
		for(j=0;j<E->nx;j++){

			// set positions for this element
			// element node numbering is anti-clockwise
			E->ele[el].node_coord[0][0] = E->wavesLon[j];
			E->ele[el].node_coord[0][1] = E->wavesLat[i]; // 0

			E->ele[el].node_coord[1][0] = E->wavesLon[j+1];
			E->ele[el].node_coord[1][1] = E->wavesLat[i]; // 1

			E->ele[el].node_coord[2][0] = E->wavesLon[j+1];
			E->ele[el].node_coord[2][1] = E->wavesLat[i+1]; // 2

			E->ele[el].node_coord[3][0] = E->wavesLon[j];
			E->ele[el].node_coord[3][1] = E->wavesLat[i+1]; // 3

			// now set the nodal value for this element
			E->ele[el].node_value[0] = E->Hs[t][i][j];
			E->ele[el].node_value[1] = E->Hs[t][i][j+1];
			E->ele[el].node_value[2] = E->Hs[t][i+1][j+1];
			E->ele[el].node_value[3] = E->Hs[t][i+1][j];

			el++;
		}
	}

    // the store_mesh call is required for the function get_owner_element
	store_mesh(E, E->nx, E->ny);
    //print_elements(E);

    // for each rho grid point
    for(i=0;i<E->nLonRho;i++){
        for(j=0;j<E->nLatRho;j++){
            // get pos of grid point
            pos[0] = E->lon_rho[i][j];
            pos[1] = E->lat_rho[i][j];

			/*
			if(pos[0] > 156.0){
				//printf("pre: i = %d, j = %d, pos = %f, %f\n",i,j,pos[0],pos[1]);
				pos[0] = 156.0;
				//printf("  modded: i = %d, j = %d, pos = %f, %f\n",i,j,pos[0],pos[1]);
			}

			if(pos[1] > -9.0){
				//printf("pre: i = %d, j = %d, pos = %f, %f\n",i,j,pos[0],pos[1]);
				pos[1] = -9.0;
				//printf("  modded: i = %d, j = %d, pos = %f, %f\n",i,j,pos[0],pos[1]);
			}
			*/

            // find out which element this lies within
			el = get_owner_element(E, pos);

            calculate_interpolation_weights(&E->ele[el], E->xi, E->eta, pos);
			interpolate_point(&E->ele[el], &E->Hs_on_roms[t][i][j]);
        }
    }
}


void interp_tp_to_roms(e *E, int t){

    int i,j;
    int el;
    double  pos[2];
    //double  interp_value;


    // setup the interpolation source grid data structures
	// source grid is the tide data
	E->nx = E->nLonWaves-1; //E->nc.x-1;
	E->ny = E->nLatWaves-1; //E->nc.y-1;
	E->nElements = E->nx*E->ny ;
	E->ele = malloc(E->nElements*sizeof(element));
	E->nodesPerEl = 4;

	E->msh.x = malloc((E->nx+1) * sizeof(double));
	E->msh.y = malloc((E->ny+1) * sizeof(double));

	init_xi_eta(E);

	el = 0;
	for(i=0;i<E->ny;i++){
		for(j=0;j<E->nx;j++){

			// set positions for this element
			// element node numbering is anti-clockwise
			E->ele[el].node_coord[0][0] = E->wavesLon[j];
			E->ele[el].node_coord[0][1] = E->wavesLat[i]; // 0

			E->ele[el].node_coord[1][0] = E->wavesLon[j+1];
			E->ele[el].node_coord[1][1] = E->wavesLat[i]; // 1

			E->ele[el].node_coord[2][0] = E->wavesLon[j+1];
			E->ele[el].node_coord[2][1] = E->wavesLat[i+1]; // 2

			E->ele[el].node_coord[3][0] = E->wavesLon[j];
			E->ele[el].node_coord[3][1] = E->wavesLat[i+1]; // 3

			// now set the nodal value for this element
			E->ele[el].node_value[0] = E->Tp[t][i][j];
			E->ele[el].node_value[1] = E->Tp[t][i][j+1];
			E->ele[el].node_value[2] = E->Tp[t][i+1][j+1];
			E->ele[el].node_value[3] = E->Tp[t][i+1][j];

			el++;
		}
	}

    // the store_mesh call is required for the function get_owner_element
	store_mesh(E, E->nx, E->ny);
    //print_elements(E);

    // for each rho grid point
    for(i=0;i<E->nLonRho;i++){
        for(j=0;j<E->nLatRho;j++){
            // get pos of grid point
            pos[0] = E->lon_rho[i][j];
            pos[1] = E->lat_rho[i][j];

			/*
			if(pos[0] > 156.0){
				//printf("pre: i = %d, j = %d, pos = %f, %f\n",i,j,pos[0],pos[1]);
				pos[0] = 156.0;
				//printf("  modded: i = %d, j = %d, pos = %f, %f\n",i,j,pos[0],pos[1]);
			}

			if(pos[1] > -9.0){
				//printf("pre: i = %d, j = %d, pos = %f, %f\n",i,j,pos[0],pos[1]);
				pos[1] = -9.0;
				//printf("  modded: i = %d, j = %d, pos = %f, %f\n",i,j,pos[0],pos[1]);
			}
			*/

            // find out which element this lies within
			el = get_owner_element(E, pos);

            calculate_interpolation_weights(&E->ele[el], E->xi, E->eta, pos);
			interpolate_point(&E->ele[el], &E->Tp_on_roms[t][i][j]);

        }
    }
}



int get_nearest_slope_index(e *E, double this_lat, double this_lon, double *lat, double *lon){

	int i;
	int the_index;
	double distance;
	double min_distance = 999999999999999.0;
	// find nearest latitude

	for(i=0;i<E->nSlopes;i++){
		distance = spheriq_dist(this_lon, this_lat, lon[i], lat[i], 0);
		if(distance < min_distance){
			min_distance = distance;
			the_index = i;
		}
	}

/*
	printf("looking for:");
	printf("\ttime = %f\n", this_time);
	printf("\tlat = %f\n", this_lat);
	printf("\tlon = %f\n", this_lon);
	printf("found this: the_index = %d\n", the_index);
	printf("\tlat = %f\n", lat[the_index]);
	printf("\tlon = %f\n", lon[the_index]);
	printf("\tsetup = %f\n", setup[0][the_index]);
	printf("\t distance = %f\n\n", min_distance);
*/
	//if(min_distance < 10000.00)
		return the_index;
	//else
	//	return 0.0;
}


double get_nearest_slope(e *E, int this_time, double this_lat, double this_lon, double **setup, double *lat, double *lon){

	int i;
	int the_index;
	double distance;
	double min_distance = 999999999999999.0;
	// find nearest latitude

	for(i=0;i<E->nSlopes;i++){
		distance = spheriq_dist(this_lon, this_lat, lon[i], lat[i], 0);
		if(distance < min_distance){
			min_distance = distance;
			the_index = i;
		}
	}

/*
	printf("looking for:");
	printf("\ttime = %f\n", this_time);
	printf("\tlat = %f\n", this_lat);
	printf("\tlon = %f\n", this_lon);
	printf("found this: the_index = %d\n", the_index);
	printf("\tlat = %f\n", lat[the_index]);
	printf("\tlon = %f\n", lon[the_index]);
	printf("\tsetup = %f\n", setup[0][the_index]);
	printf("\t distance = %f\n\n", min_distance);
*/
	//if(min_distance < 10000.00)
		return setup[this_time][the_index];
	//else
	//	return 0.0;
}

/* Linear interpolation is the simplest method of getting values at
	positions in between the data points. The points are simply joined by
	straight line segments.

	Each segment (bounded by two data points) can be interpolated independently.
	The parameter mu defines where to estimate the value on the interpolated line,
	it is 0 at the first point and 1 and the second point.

	For interpolated values between the two points mu ranges between 0 and 1.
	Values of mu outside this range result in extrapolation.
*/
double LinearInterpolate( double y1,double y2, double mu){
   return(y1*(1-mu)+y2*mu);
}

/*
Given arrays x[1..n] and y[1..n] containing a tabulated function, i.e., yi = f(xi), with
x1 < x2 <... < xN , and given values yp1 and ypn for the first derivative of the interpolating
function at points 1 and n, respectively, this routine returns an array y2[1..n] that contains
the second derivatives of the interpolating function at the tabulated points xi. If yp1 and/or
ypn are equal to 1 × 1030 or larger, the routine is signaled to set the corresponding boundary
condition for a natural spline, with zero second derivative on that boundary
*/
void spline(double *x, double *y, int n, double yp1, double ypn, double *y2) {

    double	*u;
    int		i, k;
    double	p, qn, sig, un;

	u = malloc(n*sizeof(double));

    if (yp1 > 0.99e30) {
        y2[0] = u[0] = 0.0;
    }
    else {
        y2[0] = -0.5;
        u[0] = (3.0/(x[1]-x[0]))*((y[1]-y[0])/(x[1]-x[0])-yp1);
    }

    for(i = 1; i < n-1; i++){
        sig = (x[i] - x[i-1])/(x[i+1] - x[i-1]);
        p = sig * y2[i-1] + 2.0;
        y2[i] = (sig - 1.0)/p;
        u[i] = (y[i+1] - y[i])/(x[i+1] - x[i]) - (y[i] - y[i-1])/(x[i] - x[i-1]);
        u[i] = (6.0*u[i]/(x[i+1] - x[i-1]) - sig*u[i-1])/p;

    }

    if (ypn > 0.99e30) {
        qn = un = 0.0;
    }
    else {
        qn = 0.5;
        un = (3.0/(x[n-1] - x[n-2]))*(ypn - (y[n-1] - y[n-2])/(x[n-1] - x[n-2]));
    }

    y2[n-1] = (un - qn*u[n-2])/(qn*y2[n-2] + 1.0);

    for (k = n-2; k >= 0; k--){
        y2[k] = y2[k] * y2[k+1] + u[k];
    }

	free(u);

}


double splint(double *xa, double *ya, double *y2a, int n, double x) {
	int klo,khi,k;
	double h,b,a;
	static int pklo=0,pkhi=1;

	if (xa[pklo] <= x && xa[pkhi] > x) {
		klo = pklo;
		khi = pkhi;
	} else {
		klo = 0;
		khi = n - 1;
		while (khi - klo > 1) {
			k = (khi + klo) >> 1;
			if (xa[k] > x) {
				khi = k;
			} else {
				klo = k;
			}
		}
	}

	h = xa[khi] - xa[klo];
	if (h == 0) {
		fprintf(stderr,"Bad xa input to function splint()\n");
		exit(1);
	}

	a = (xa[khi] - x)/h;
	b = (x - xa[klo])/h;
	return a*ya[klo] + b*ya[khi] + ((a*a*a - a)*y2a[klo] + (b*b*b - b)*y2a[khi])*(h*h)/6.0;
}

/* send this function the control points and an array of desired interpolated
	points within the control points and it will return the interpolated
	values at those points */
void time_interp_field( 	double *xpts, double *ypts, int npts_in,
								double *interp_x, double * interp_y, int interp_npts){
	int		i;
	double	*interp_spline;
	// setup the spline stuff
	interp_spline = malloc(npts_in*sizeof(double));
	//void spline(double *x, double *y, int n, double yp1, double ypn, double *y2)
	// x is the x axis data
	// y is the value at each x
	// npts is how many control points we have
	// yp1 is the derivative at x = 0
	// ypn is the derivative at x = n
	// *y2 is the spline
	spline(&xpts[0], &ypts[0], npts_in, 0.0, 0.0, &interp_spline[0]);

	//for(i=0;i<npts_in;i++)
	//	printf("xpts[%d] = %f,  =ypts[%d] = %f\n",i,xpts[i],i,ypts[i]);

	// now do the spline interpolation for this time
	for(i=0;i<interp_npts;i++){
		/*
		Given the arrays xa[1..n] and ya[1..n], which tabulate a function (with the xai’s in order),
		and given the array y2a[1..n], which is the output from spline above, and given a value of
		x, this routine returns a cubic-spline interpolated value y
		*/
		//printf("interp_x[%d] = %f\n", i, interp_x[i]);
		interp_y[i] = splint(&xpts[0], &ypts[0], &interp_spline[0], npts_in, interp_x[i]);
		//printf("interp value is %f\n", interp_y[i]);

	}

	free(interp_spline);

}
