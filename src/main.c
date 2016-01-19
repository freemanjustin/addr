// gridr
//
// freeman.justin@gmail.com


#include "grid.h"

int main(int argc,char **argv)
{
	e	*E;
    int     Lm,Mm,Lp,Mp, L, M;
    int     i,j;
    double  latdist;
    double  height, length;

    double  el, xl;


	// addr netcdf variables
	int	ncid;
	int varid;
	int retval;
	size_t attlen = 0;

	// nn variables
	int	n_valid_points;

	// for time conversion
	static	char *calendar = "Standard";
    ut_system	*u_system;
    ut_unit	*ref_time;
    double	tval_inv, sec;
    int	ierr, yr, mo, day, hr, min;

	// malloc the work struct
	E = malloc(sizeof(e));
	if(E==NULL) fail("Malloc failed\n");

	// parse command line arguments
	if(argc < 3){
		fail("need an input xml file and an output filename\n");
	}
	else{
		get_command_line_arg_as_string(&E->input_xml, argv[1]);
		get_command_line_arg_as_string(&E->roms_input, argv[2]);
		get_command_line_arg_as_string(&E->tide_input, argv[3]);
		get_command_line_arg_as_string(&E->fname, argv[4]);
	}

	// initialize the time converison libs
	// Initialize the udunits-2 library
    ut_set_error_message_handler(ut_ignore);
    if( (u_system = ut_read_xml( NULL )) == NULL ) {
        fprintf( stderr, "Error initializing udunits-2 unit system\n" );
        exit(-1);
        }
    ut_set_error_message_handler(ut_write_to_stderr);

    // read input XML
    //get_params(E);

	printf("roms file = %s\n", E->roms_input);
	printf("tide file = %s\n", E->tide_input);
	printf("outputfile = %s\n", E->fname);

	// read the target grid - this is the roms output curvilinear grid
	// in this code we are interpolating the tide data onto the rho points

	/* roms output data required

	double lon_rho(eta_rho, xi_rho) ;
		lon_rho:long_name = "longitude of RHO-points" ;
		lon_rho:units = "degree_east" ;
		lon_rho:standard_name = "longitude" ;
		lon_rho:field = "lon_rho, scalar" ;
	double lat_rho(eta_rho, xi_rho) ;
		lat_rho:long_name = "latitude of RHO-points" ;
		lat_rho:units = "degree_north" ;
		lat_rho:standard_name = "latitude" ;
		lat_rho:field = "lat_rho, scalar" ;
	*/
	// open the roms_input file
	// open the file
    if((retval = nc_open(E->roms_input, NC_NOWRITE, &ncid)))
        fail("failed to open roms input file: error is %d\n",retval);

	// get the time data
	if((retval = nc_inq_dimid(ncid, "ocean_time", &varid)))
        fail("failed to get roms ocean_time dimid: error is %d\n",retval);

    if((retval = nc_inq_dimlen(ncid,varid,&E->nTimeRoms)))
        fail("failed to get roms lat dimid: error is %d\n",retval);

    printf("ocean_time = %zu\n", E->nTimeRoms);

    // get the lat dimension sizes
    if((retval = nc_inq_dimid(ncid, "xi_rho", &varid)))
        fail("failed to get roms lat dimid: error is %d\n",retval);

    if((retval = nc_inq_dimlen(ncid,varid,&E->nLatRho)))
        fail("failed to get roms lat dimid: error is %d\n",retval);

    printf("xi_rho = %zu\n", E->nLatRho);

    // get the lon dimension sizes
    if((retval = nc_inq_dimid(ncid, "eta_rho", &varid)))
        fail("failed to get roms lon_rho dimid: error is %d\n",retval);

    if((retval = nc_inq_dimlen(ncid,varid,&E->nLonRho)))
        fail("failed to read roms lon_rho data: error is %d\n",retval);

    //cout << "nlats = " << nlat << endl;
    //cout << "nlons = " << nlon << endl;

    printf("eta_rho = %zu\n", E->nLonRho);

    // malloc room for the arrays
	E->romsTime = malloc(E->nTimeRoms*sizeof(double));
	E->lat_rho = malloc2d_double(E->nLatRho, E->nLonRho);
	E->lon_rho = malloc2d_double(E->nLatRho, E->nLonRho);
	E->mask_rho = malloc2d_double(E->nLatRho, E->nLonRho);


    // read the data
	nc_inq_varid(ncid, "ocean_time", &varid);
    if((retval = nc_get_var_double(ncid, varid, &E->romsTime[0])))
		fail("failed to read roms ocean_time data: error is %d\n", retval);

	// get the time metadata units
	nc_inq_attlen (ncid, varid, "units", &attlen);
	E->roms_time_units = (char *) malloc(attlen + 1);  /* + 1 for trailing null */
	nc_get_att_text(ncid, varid, "units", E->roms_time_units);
	printf("units = %s\n", E->roms_time_units);

	printf("romsTime[0] = %f\n", E->romsTime[0]);
	// Make the Calendar calls
    //tval = 86460.0;	/* in seconds, this is 1 day and 1 minute */
    //tval = 8580;

	/* Parse the units strings */
    if( (ref_time = ut_parse( u_system, E->roms_time_units, UT_ASCII )) == NULL ) {
        fprintf( stderr, "Error parsing units string \"%s\"\n", E->roms_time_units );
        exit(-1);
        }

    if( (ierr = utCalendar2_cal( E->romsTime[0], ref_time, &yr, &mo, &day, &hr, &min, &sec, calendar )) != 0 ) {
        fprintf( stderr, "Error on utCalendar2_cal call: %d\n", ierr );
        exit(-1);
        }
    printf( "this date is %04d-%02d-%02d %02d:%02d:%06.3lf\n",yr, mo, day, hr, min, sec );


    nc_inq_varid(ncid, "lat_rho", &varid);
    if((retval = nc_get_var_double(ncid, varid, &E->lat_rho[0][0])))
		fail("failed to read roms lat_rho data: error is %d\n", retval);

		printf("lat_rho[0][0] = %f\n", E->lat_rho[0][0]);

    nc_inq_varid(ncid, "lon_rho", &varid);
    if((retval = nc_get_var_double(ncid, varid, &E->lon_rho[0][0])))
		fail("failed to read roms lat_rho data: error is %d\n", retval);

	printf("lon_rho[0][0] = %f\n", E->lon_rho[0][0]);

	// get the rho_mask

	nc_inq_varid(ncid, "mask_rho", &varid);
    if((retval = nc_get_var_double(ncid, varid, &E->mask_rho[0][0])))
		fail("failed to read roms mask_rho data: error is %d\n", retval);

	// close the roms file for now
	// will probably open the roms file as read-write so we can dump
	// the interpolated tides onto it
	nc_close(ncid);


	// read in the tide file

	if((retval = nc_open(E->tide_input, NC_NOWRITE, &ncid)))
        fail("failed to open tide input file: error is %d\n",retval);

    // get the lat dimension sizes
    if((retval = nc_inq_dimid(ncid, "lat", &varid)))
        fail("failed to get tide lat dimid: error is %d\n",retval);

    if((retval = nc_inq_dimlen(ncid,varid,&E->nLatTide)))
        fail("failed to get tide lat dim length: error is %d\n",retval);

    printf("lat = %zu\n", E->nLatTide);

    // get the lon dimension sizes
    if((retval = nc_inq_dimid(ncid, "lon", &varid)))
        fail("failed to get tide lon dimid: error is %d\n",retval);

    if((retval = nc_inq_dimlen(ncid,varid,&E->nLonTide)))
        fail("failed to read tide lon dim length: error is %d\n",retval);

    printf("tide lon = %zu\n", E->nLonTide);

	// get the lev dimension sizes
    if((retval = nc_inq_dimid(ncid, "lev", &varid)))
        fail("failed to get tide lev dimid: error is %d\n",retval);

    if((retval = nc_inq_dimlen(ncid,varid,&E->nLevTide)))
        fail("failed to read tide lon data: error is %d\n",retval);

    printf("tide lev = %zu\n", E->nLevTide);

	// get the time dimension sizes
    if((retval = nc_inq_dimid(ncid, "time", &varid)))
        fail("failed to get tide tide dimid: error is %d\n",retval);

    if((retval = nc_inq_dimlen(ncid,varid,&E->nTimeTide)))
        fail("failed to read tide time data: error is %d\n",retval);

    printf("tide time = %zu\n", E->nTimeTide);

    // malloc room for the arrays
	//(time, lev, lat, lon)
	E->tideTime = malloc(E->nTimeTide*sizeof(double));
	E->tideLon = malloc(E->nLonTide*sizeof(double));
	E->tideLat = malloc(E->nLatTide*sizeof(double));
    E->tide_data = malloc4d_double(E->nTimeTide, E->nLevTide, E->nLatTide, E->nLonTide);


	// read in the time data
	nc_inq_varid(ncid, "time", &varid);
    if((retval = nc_get_var_double(ncid, varid, &E->tideTime[0])))
		fail("failed to read tide time data data: error is %d\n",retval);
	// get the time metadata units
	nc_inq_attlen (ncid, varid, "units", &attlen);
	E->tide_time_units = (char *) malloc(attlen + 1);  /* + 1 for trailing null */
	nc_get_att_text(ncid, varid, "units", E->tide_time_units);
	lr_pack(E->tide_time_units);
	printf("tide time units = %s\n", E->tide_time_units);

	printf("tideTime[0] = %f\n", E->tideTime[0]);
	// Make the Calendar calls
	/* Parse the units strings */
    if( (ref_time = ut_parse( u_system, E->tide_time_units, UT_ASCII )) == NULL ) {
        fprintf( stderr, "Error parsing units string \"%s\"\n", E->roms_time_units );
        exit(-1);
        }

    if( (ierr = utCalendar2_cal( E->tideTime[0], ref_time, &yr, &mo, &day, &hr, &min, &sec, calendar )) != 0 ) {
        fprintf( stderr, "Error on utCalendar2_cal call: %d\n", ierr );
        exit(-1);
        }
    printf( "this date is %04d-%02d-%02d %02d:%02d:%06.3lf\n",yr, mo, day, hr, min, sec );





    // read the data
	nc_inq_varid(ncid, "lat", &varid);
    if((retval = nc_get_var_double(ncid, varid, &E->tideLat[0])))
		fail("failed to read tide lat data data: error is %d\n",retval);

	nc_inq_varid(ncid, "lon", &varid);
    if((retval = nc_get_var_double(ncid, varid, &E->tideLon[0])))
		fail("failed to read tide lon data data: error is %d\n",retval);

    nc_inq_varid(ncid, "z", &varid);
    if((retval = nc_get_var_double(ncid, varid, &E->tide_data[0][0][0][0])))
		fail("failed to read tide data data: error is %d\n",retval);

	// close the tide file
	nc_close(ncid);

		// do a natural neighbour interpolation on the tide data we just read in
		// to fill in the land masked values before interpolating onto the ROMS grid
		// actually, we can probably skip the bilinear inperpolation step and
		// natural neighbor interpolate directly onto the ROMS grid. Since this is a
		// test run, then do it in two steps...

		// find out how many valid data points we have
		// and setup the input array for lib-nn
		E->nn_diff = malloc(E->nLonTide * E->nLatTide * sizeof(point));
		E->nn_n = 0;
		for(i=0;i<E->nLatTide;i++){
			for(j=0;j<E->nLonTide;j++){
				if(E->tide_data[0][0][i][j] > -999.0){
					E->nn_diff[E->nn_n].y = E->tideLat[i];
					E->nn_diff[E->nn_n].x = E->tideLon[j];
					E->nn_diff[E->nn_n].z = E->tide_data[0][0][i][j];

					E->nn_n++;
				}
			}
		}

		// now do the interpolation via lib-nn
		E->nn_weight = 0.0;
		// so have E->weight = 0.0; somewhere
		// E->diff is the input data list
		// E->interp is the output data list
		// E-> nx and E->ny is the grid dimensions
		E->nn_dx = fabs(E->tideLat[1] - E->tideLat[0]); // lon
		E->nn_dy = fabs(E->tideLon[1] - E->tideLon[0]); // lat
		// do the natural neighbour interpolation
		// this will interpolate our scattered point data
		// onto an orthogonal grid with dimensions
		// E->nx x E->ny
		// figure out E->nx and E->ny
		get_mesh_dimensions(E, E->nn_n, E->nn_diff);
		// and malloc room for the interpolated data
		// we already know the mesh dimensions from the input tide file
		// E->nx = number of longitudes in file
		// E->ny = number of latitudes in file
		E->nn_interp = malloc(E->nLonTide * E->nLatTide * sizeof(point));
		printf("doing nn interpolation...\n"); fflush(stdout);
		// E->n is the number of points in E->diff
		// E->weight is a constant (equal to 0.0 and declared in the input xml file)


		nn_interp_to_mesh(E, E->nn_n, E->nn_weight, E->nn_diff, E->nLonTide, E->nLatTide, E->nn_interp);
		printf("done\n"); fflush(stdout);

		// put the interpolated data into the tide array so we can write it out to
		// a netcdf file
		double *nn_interp;
		nn_interp = malloc(E->nn_nx * E->nn_ny * sizeof(double));
		printf("nn interp field: nn_nx = %d, nn_ny = %d\n", E->nn_nx, E->nn_ny);
		for(i=0;i<E->nn_nx*E->nn_ny;i++){
			//fprintf(fptr,"%.15g %.15g %.15g\n", p[i].x, p[i].y, p[i].z);
			nn_interp[i] = E->nn_interp[i].z;
		}

	// write out the interped field
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

	// interpolate the tide data for each lon_rho and lat_rho point

	// target grid for tide data
	// this needs a time dimension from the roms file
	E->tide_on_roms = malloc2d_double(E->nLatRho, E->nLonRho);
	for(i=0;i<E->nLatRho;i++)
		for(j=0;j<E->nLonRho;j++)
			E->tide_on_roms[i][j] = 0.0;


	// temporarily splat the nn_interped field over the original data
	int count = 0;
	for(i=0;i<E->nLatTide;i++){
		for(j=0;j<E->nLonTide;j++){
			E->tide_data[0][0][i][j] = nn_interp[count];
			count++;
		}
	}

	// interpolate tide to roms grid
    interp_tide_to_roms(E);


	// apply land-sea mask
	for(i=0;i<E->nLatRho;i++){
		for(j=0;j<E->nLonRho;j++){
			if(E->mask_rho[i][j] == 0)
				E->tide_on_roms[i][j] = NC_FILL_DOUBLE;
		}
	}
	// write the interpolated field to file
	write_netcdf(E);


	return 0;
}
