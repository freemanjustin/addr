// gridr
//
// freeman.justin@gmail.com



#include "grid.h"



void write_netcdf(e *E){

	//float	fillValue = -1e34;
	float fillValue = NC_FILL_FLOAT;


	create_netcdf(E, E->fname, &E->ncid);

    defdims_netcdf(E);

	// setup variables
	defvars(E);


	/*
	// add global metadata
	add_global_metadata(E, E->ncid);
	*/


    // write the data to the netcdf file
    write_data(E);

    nc_close(E->ncid);
}


void create_netcdf(e *E, char *fname, int *ncid){

	//int	old_fill_mode;

	if ( (E->retval = nc_create(fname, NC_CLOBBER, ncid) ) )
		fail("couldn't create netcdf out file %s\n",fname);

	// set the fill mode to be NO_FILL
	//nc_set_fill(*ncid, NC_NOFILL, &old_fill_mode);

}

void defdims_netcdf(e *E){

	E->one = 1;

	// define the dimensions
	if ((E->retval = nc_def_dim(E->ncid, "xi_rho", E->nLatRho, &E->xi_rho_dimid)))
		fail("nc_def_dim failed!\n");

    if ((E->retval = nc_def_dim(E->ncid, "eta_rho", E->nLonRho, &E->eta_rho_dimid)))
		fail("nc_def_dim failed!\n");


	// end define mode for this file
	if ((E->retval = nc_enddef(E->ncid)))
      	fail("nc_enddef failed\n");

}


void defvars(e *E){

    // setup dimids

    E->dimIdsRho[0] = E->eta_rho_dimid;
    E->dimIdsRho[1] = E->xi_rho_dimid;


	defvar_netcdf(E, E->ncid, "tide", NC_DOUBLE, 2, &E->dimIdsRho[0], &E->vid_tide);
    add_txt_attribute_netcdf(E, E->ncid, E->vid_h, "long_name", "tide at RHO-points");
	add_txt_attribute_netcdf(E, E->ncid, E->vid_h, "units", "m");

	defvar_netcdf(E, E->ncid, "lat_rho", NC_DOUBLE, 2, E->dimIdsRho, &E->vid_lat_rho);
    add_txt_attribute_netcdf(E, E->ncid, E->vid_lat_rho, "long_name", "latitude at RHO-points");
	add_txt_attribute_netcdf(E, E->ncid, E->vid_lat_rho, "units", "degree_north");
    add_txt_attribute_netcdf(E, E->ncid, E->vid_lat_rho, "_CoordinateAxisType", "Lat");

	defvar_netcdf(E, E->ncid, "lon_rho", NC_DOUBLE, 2, &E->dimIdsRho[0], &E->vid_lon_rho);
    add_txt_attribute_netcdf(E, E->ncid, E->vid_lon_rho, "long_name", "longitude at RHO-points");
	add_txt_attribute_netcdf(E, E->ncid, E->vid_lon_rho, "units", "degree_east");
    add_txt_attribute_netcdf(E, E->ncid, E->vid_lon_rho, "_CoordinateAxisType", "Lon");

}



void defvar_netcdf(e *E, int ncid, char *var_name, nc_type type, int dims, int *dimIds, int *varid){


	if((E->retval = nc_redef(ncid)))
		fail("nc_redef failed\n");

	if ((E->retval = nc_def_var(ncid, var_name, type, dims, dimIds, varid)))
		fail("nc_def_var failed\n");

	if ((E->retval = nc_enddef(ncid)))
		fail("nc_enddef failed\n");

}


void add_txt_attribute_netcdf(e *E, int ncid, int varid, char* att_name, char* att_value){

	if((E->retval = nc_redef(ncid)))
		fail("nc_redef failed\n");

	if(( E->retval = nc_put_att_text(ncid, varid, att_name, strlen(att_value), att_value)))
		fail("nc_put_att_txt failed\n");

	if ((E->retval = nc_enddef(ncid)))
		fail("nc_enddef failed\n");
}


void add_global_metadata(e *E, int ncid){

    // TODO!
}


void write_data(e *E){


    // write tide
    if ((E->retval = nc_put_var_double(E->ncid, E->vid_tide, &E->tide_on_roms[0][0])))
        fail("put_var_ failed. Error code = %d\n",E->retval);


    // lat_rho
    if ((E->retval = nc_put_var_double(E->ncid, E->vid_lat_rho, &E->lat_rho[0][0])))
        fail("put_var_ failed for lat_rho. Error code = %d\n",E->retval);

    // lon_rho
    if ((E->retval = nc_put_var_double(E->ncid, E->vid_lon_rho, &E->lon_rho[0][0])))
        fail("put_var_ failed. Error code = %d\n",E->retval);


}
