// gridr
//
// freeman.justin@gmail.com


#include "grid.h"

int main(int argc,char **argv)
{
	e	*E;
    int     i,j,t;
	int count;


	// addr netcdf variables
	int	ncid;
	int varid;
	int retval;
	size_t attlen = 0;

	// malloc the work struct
	E = malloc(sizeof(e));
	if(E==NULL) fail("Malloc failed\n");

	// parse command line arguments
	if(argc < 3){
		printf("usage:\n");
		printf("\taddr [roms history netcdf file] [tide netcdf file] [auswave netcdf file] [output filename]\n");
		printf("\n\tauswave data must include Hs and Tp\n");
		exit(1);
	}
	else{
		get_command_line_arg_as_string(&E->roms_input, argv[1]);
		get_command_line_arg_as_string(&E->tide_input, argv[2]);
		get_command_line_arg_as_string(&E->wave_input, argv[3]);
		get_command_line_arg_as_string(&E->fname, argv[4]);
	}

	//printf("roms file = %s\n", E->roms_input);
	//printf("tide file = %s\n", E->tide_input);
	//printf("wave setup file = %s\n", E->wave_input);
	//printf("outputfile = %s\n", E->fname);

	// initialize the time converison libs
	// Initialize the udunits-2 library

	sprintf(E->calendar,"Standard");
    ut_set_error_message_handler(ut_ignore);
    if( (E->u_system = ut_read_xml( NULL )) == NULL ) {
        fprintf( stderr, "Error initializing udunits-2 unit system\n" );
        exit(-1);
        }
    ut_set_error_message_handler(ut_write_to_stderr);

	// read in the roms input data
	process_roms(E);

	// get setup on roms grid
	process_auswave(E);

	// get tide on roms grid
	process_tides(E);

	// add the components
	add(E);

	// write the interpolated field to file
	//write_netcdf(E);
	//printf("writing coastal data...");
	write_coastal_data(E);
	//printf("done\n");

	// write out time series data for each coastal point
	#ifdef CHECK
	write_time_series(E);
	#endif

	/*
	// write out the interped field coming from the natural neighbor interpolation
	//E->nn_nx * E->nn_ny, E->nn_interp,
	int lat_dimid, lon_dimid, dimIds[2];
	int lat_varid, lon_varid, tide_interp_varid;
	// create the file
	nc_create("tide_interp.nc", NC_CLOBBER, &ncid);
	// def dimensions
	nc_def_dim(ncid, "lat", E->nn_ny, &lat_dimid);
	nc_def_dim(ncid, "lon", E->nn_nx, &lon_dimid);
	// def vars
	dimIds[0] = lat_dimid;
	dimIds[1] = lon_dimid;
	nc_def_var(ncid, "lat", NC_DOUBLE, 1, &dimIds[0], &lat_varid);
	nc_def_var(ncid, "lon", NC_DOUBLE, 1, &dimIds[1], &lon_varid);
	nc_def_var(ncid, "tide_interp", NC_DOUBLE, 2, dimIds, &tide_interp_varid);
	nc_enddef(ncid);
	// write the data
	nc_put_var_double(ncid, lat_varid, &E->tideLat[0]);
	nc_put_var_double(ncid, lon_varid, &E->tideLon[0]);
	nc_put_var_double(ncid, tide_interp_varid, nn_interp);
	// close the file
	nc_close(ncid);
	*/


	return 0;
}
