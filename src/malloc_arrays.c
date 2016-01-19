// gridr
//
// freeman.justin@gmail.com


#include "grid.h"

void malloc_arrays(e *E){

    E->lat_rho = malloc2d_double(E->g.nY+1, E->g.nX);
    E->lon_rho = malloc2d_double(E->g.nY+1, E->g.nX);

}
