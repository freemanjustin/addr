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

	// nn variables
	double *nn_interp;

	// for time conversion
	static	char *calendar = "Standard";
    ut_system	*u_system;
    double	sec;
    int	ierr, yr, mo, day, hr, min;

	// malloc the work struct
	E = malloc(sizeof(e));
	if(E==NULL) fail("Malloc failed\n");

	// parse command line arguments
	if(argc < 3){
		fail("need an input xml file and an output filename\n");
	}
	else{
		get_command_line_arg_as_string(&E->roms_input, argv[1]);
		get_command_line_arg_as_string(&E->tide_input, argv[2]);
		get_command_line_arg_as_string(&E->wave_input, argv[3]);
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

	printf("roms file = %s\n", E->roms_input);
	printf("tide file = %s\n", E->tide_input);
	printf("wave setup file = %s\n", E->wave_input);
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

    printf("eta_rho = %zu\n", E->nLonRho);

    // malloc room for the arrays
	E->romsTime = malloc(E->nTimeRoms*sizeof(double));
	E->lat_rho = malloc2d_double(E->nLonRho, E->nLatRho);
	E->lon_rho = malloc2d_double(E->nLonRho, E->nLatRho);
	E->mask_rho = malloc2d_double(E->nLonRho, E->nLatRho);
	E->zeta = malloc3d_double(E->nTimeRoms, E->nLonRho, E->nLatRho);


    // read the data from the ROMS output file
	nc_inq_varid(ncid, "ocean_time", &varid);
    if((retval = nc_get_var_double(ncid, varid, &E->romsTime[0])))
		fail("failed to read roms ocean_time data: error is %d\n", retval);



	// normalize the time information between the roms_his file
	// and the tide data file
	// get everything in ROMS ocean_time
	// get the time metadata units
	nc_inq_attlen (ncid, varid, "units", &attlen);
	E->roms_time_units = (char *) malloc(attlen + 1);  /* + 1 for trailing null */
	E->roms_time_units[attlen] = '\x0';
	nc_get_att_text(ncid, varid, "units", E->roms_time_units);
	printf("units = %s\n", E->roms_time_units);

	printf("romsTime[0] = %f\n", E->romsTime[0]);
	// Make the Calendar calls
    //tval = 86460.0;	/* in seconds, this is 1 day and 1 minute */
    //tval = 8580;

	/* Parse the units strings */
    if( (E->roms_ref_time = ut_parse( u_system, E->roms_time_units, UT_ASCII )) == NULL ) {
        fprintf( stderr, "Error parsing units string \"%s\"\n", E->roms_time_units );
        exit(-1);
        }

    if( (ierr = utCalendar2_cal( E->romsTime[0], E->roms_ref_time, &yr, &mo, &day, &hr, &min, &sec, calendar )) != 0 ) {
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

	// get the sea_surface_height
	nc_inq_varid(ncid, "zeta", &varid);
    if((retval = nc_get_var_double(ncid, varid, &E->zeta[0][0][0])))
		fail("failed to read roms height data: error is %d\n", retval);

	// close the roms file for now
	// will probably open the roms file as read-write so we can dump
	// the interpolated tides onto it
	nc_close(ncid);


















	// read in the wave setup data
	// open the file
    if((retval = nc_open(E->wave_input, NC_NOWRITE, &ncid)))
        fail("failed to open wave setup input file: error is %d\n",retval);

	// get the time data
	if((retval = nc_inq_dimid(ncid, "time", &varid)))
        fail("failed to get wave setup dimid: error is %d\n",retval);

    if((retval = nc_inq_dimlen(ncid,varid,&E->nTimeWaves)))
        fail("failed to get roms lat dimlen: error is %d\n",retval);

    printf("waves_times = %zu\n", E->nTimeWaves);

    // get the lat dimension sizes
    if((retval = nc_inq_dimid(ncid, "station", &varid)))
        fail("failed to get waves station dimid: error is %d\n",retval);

    if((retval = nc_inq_dimlen(ncid,varid,&E->nStationWaves)))
        fail("failed to get waves statuib dimlen: error is %d\n",retval);

    printf("waves station = %zu\n", E->nStationWaves);

    // malloc room for the arrays
	E->wavesLat = malloc(E->nStationWaves*sizeof(double));
	E->wavesLon = malloc(E->nStationWaves*sizeof(double));
	E->wavesTime = malloc(E->nTimeWaves*sizeof(double));
	E->setup = malloc2d_double(E->nTimeWaves, E->nStationWaves);


    // read the data from the waves output file
	nc_inq_varid(ncid, "time", &varid);
    if((retval = nc_get_var_double(ncid, varid, &E->wavesTime[0])))
		fail("failed to read waves time data: error is %d\n", retval);

	// normalize the time information between the roms_his file
	// and the tide data file
	// get everything in ROMS ocean_time
	// get the time metadata units
	nc_inq_attlen (ncid, varid, "units", &attlen);
	E->waves_time_units = (char *) malloc(attlen + 1);  /* + 1 for trailing null */
	E->waves_time_units[attlen] = '\x0';
	nc_get_att_text(ncid, varid, "units", E->waves_time_units);
	printf("waves time units = %s\n", E->waves_time_units);

	printf("wavesTime[0] = %f\n", E->wavesTime[0]);
	// Make the Calendar calls
    //tval = 86460.0;	/* in seconds, this is 1 day and 1 minute */
    //tval = 8580;


	// Parse the units strings
    if( (E->waves_ref_time = ut_parse( u_system, E->waves_time_units, UT_ASCII )) == NULL ) {
        fprintf( stderr, "Error parsing units string \"%s\"\n", E->waves_time_units );
        exit(-1);
        }

	/*
    if( (ierr = utCalendar2_cal( E->wavesTime[0], E->waves_ref_time, &yr, &mo, &day, &hr, &min, &sec, calendar )) != 0 ) {
        fprintf( stderr, "Error on utCalendar2_cal call: %d\n", ierr );
        exit(-1);
        }
    printf( "this date is %04d-%02d-%02d %02d:%02d:%06.3lf\n",yr, mo, day, hr, min, sec );
	*/

	// put the tide time on the same time units as the roms time
	for(t=0;t<E->nTimeWaves;t++){

		//printf("PRE: waveTIme[t] = %f  ", E->wavesTime[t]);
		// convert tide time into year, month, day, hour, minute and seconds
	    if( (ierr = utCalendar2_cal( E->wavesTime[t], E->waves_ref_time, &yr, &mo, &day, &hr, &min, &sec, calendar )) != 0 ) {
	        fprintf( stderr, "Error on utCalendar2_cal call: %d\n", ierr );
	        exit(-1);
	        }
	    //printf( "this date is %04d-%02d-%02d %02d:%02d:%06.3lf\n",yr, mo, day, hr, min, sec );

		// convert this date to be on the same units as the roms time
		if( (ierr = utInvCalendar2_cal( yr, mo, day, hr, min, sec, E->roms_ref_time, &E->wavesTime[t], calendar )) != 0 ) {
			fprintf( stderr, "Error on utCalendar2_cal call: %d\n", ierr );
			exit(-1);
			}
		//printf( "POST: %04d-%02d-%02d %02d:%02d:%06.3lf is %lf %s in the %s calendar\n",
		//	yr, mo, day, hr, min, sec, E->wavesTime[t], E->roms_time_units, calendar );
	}



	nc_inq_varid(ncid, "lat", &varid);
    if((retval = nc_get_var_double(ncid, varid, &E->wavesLat[0])))
		fail("failed to read waves lat data: error is %d\n", retval);

		printf("waves lat[0] = %f\n", E->wavesLat[0]);

    nc_inq_varid(ncid, "lon", &varid);
    if((retval = nc_get_var_double(ncid, varid, &E->wavesLon[0])))
		fail("failed to read waves lon data: error is %d\n", retval);

	printf("waves lon[0] = %f\n", E->wavesLon[0]);

	// get the rho_mask
	nc_inq_varid(ncid, "setup", &varid);
    if((retval = nc_get_var_double(ncid, varid, &E->setup[0][0])))
		fail("failed to read waves setup data: error is %d\n", retval);

	printf("waves setup[0][0] = %f\n", E->setup[0][0]);

	// close the file
	nc_close(ncid);


	// malloc room for the output nearest neighbor interp wave setup data
	// target grid for tide data

	// jNOTE: fix up the size of the time dimension here!
	E->setup_on_roms = malloc3d_double(E->nTimeWaves, E->nLonRho, E->nLatRho);
	// malloc room for the time interpolated data
	E->setup_on_roms_time_interp = malloc3d_double(E->nTimeRoms, E->nLonRho, E->nLatRho);

	E->coastline_mask = malloc2d_double(E->nLonRho, E->nLatRho);
	double **padded_mask = malloc2d_double(E->nLonRho+2, E->nLatRho+2);
	// initialize the field to be fill_value everywhere
	printf("i size nLatRho = %d\n", E->nLatRho);
	printf("j size nLonRho = %d\n", E->nLonRho);

	for(i=0;i<E->nLonRho+2;i++){
		for(j=0;j<E->nLatRho+2;j++){
			//E->setup_on_roms[i][j] = -1;//E->mask_rho[i][j];
			if(i==0)
				padded_mask[i][j] = 1.0;
			if(j==0)
				padded_mask[i][j] = 1.0;
			if(i==E->nLonRho+1)
				padded_mask[i][j] = 1.0;
			if(j==E->nLatRho+1)
				padded_mask[i][j] = 1.0;
		}
	}

	for(i=0;i<E->nLonRho;i++){
		for(j=0;j<E->nLatRho;j++){
			padded_mask[i+1][j+1] = E->mask_rho[i][j];
		}
	}


	double sum;
	for(i=1;i<E->nLonRho-1;i++){
		for(j=1;j<E->nLatRho-1;j++){
				sum = 0.0;
				sum += padded_mask[i-1][j-1] * 1.0/8.0;
				sum += padded_mask[i-1][j] * 1.0/8.0;
				sum += padded_mask[i-1][j+1] * 1.0/8.0;

				sum += padded_mask[i][j-1] * 1.0/8.0;
				sum += padded_mask[i][j] * 0;
				sum += padded_mask[i][j+1] * 1.0/8.0;

				sum += padded_mask[i+1][j-1] * 1.0/8.0;
				sum += padded_mask[i+1][j] * 1.0/8.0;
				sum += padded_mask[i+1][j+1] * 1.0/8.0;

				//printf("i = %d, j = %d, sum = %f\n",i,j,sum);
				if (sum == 1)	// all ocean neighbors
					E->coastline_mask[i-1][j-1] = 0;
				else if(sum < 1)
					E->coastline_mask[i-1][j-1] = 1;	// is a coast

				E->coastline_mask[i-1][j-1] = E->coastline_mask[i-1][j-1] && E->mask_rho[i-1][j-1];
		}
	}

	free(padded_mask);

	// get coastal zeta from roms_his
	double	***zeta_coast = malloc3d_double(E->nTimeRoms, E->nLonRho, E->nLatRho);
	for(t=0;t<E->nTimeRoms;t++){
		for(i=0;i<E->nLonRho;i++){
			for(j=0;j<E->nLatRho;j++){
				zeta_coast[t][i][j] = E->zeta[t][i][j] * E->coastline_mask[i][j];
			}
		}
	}

	// assign closest wavesetup value to costline derived from the roms rho_mask
	// E->setup = malloc2d_double(E->nTimeWaves, E->nStationWaves);
	double this_time;
	double this_lat;
	double this_lon;
	int	**nearest_index = malloc2d_int(E->nLonRho, E->nLatRho);

	// first get the index mapping for each coastal cell
	for(i=0;i<E->nLonRho;i++){
		for(j=0;j<E->nLatRho;j++){
			if(E->coastline_mask[i][j] == 1){ // if this point is a coastal point
				// get longitude and latitude of the point
				this_time = t;
				this_lat = E->lat_rho[i][j];
				this_lon = E->lon_rho[i][j];
				nearest_index[i][j] = get_nearest_setup_index(E, this_time, this_lat, this_lon, E->setup, E->wavesLat, E->wavesLon);
				//exit(1);

			}
			else{
				//printf("fill it: i = %d, j = %d\n",i,j);
				nearest_index[i][j] = -999;
			}
		}
	}


	for(t=0;t<E->nTimeWaves;t++){
		printf("#### t = %d\n",t);
		for(i=0;i<E->nLonRho;i++){
			for(j=0;j<E->nLatRho;j++){
				if(E->coastline_mask[i][j] == 1){ // if this point is a coastal point
					/*
					// get longitude and latitude of the point
					this_time = t;
					this_lat = E->lat_rho[i][j];
					this_lon = E->lon_rho[i][j];
					E->setup_on_roms[t][i][j] = get_nearest_setup(E, this_time, this_lat, this_lon, E->setup, E->wavesLat, E->wavesLon);
					//exit(1);
					*/
					E->setup_on_roms[t][i][j] = E->setup[t][nearest_index[i][j]];
				}
				else{
					//printf("fill it: i = %d, j = %d\n",i,j);
					E->setup_on_roms[t][i][j] = 0.0;//NC_FILL_DOUBLE;
				}
			}
		}
	}

	// time interpolate the wavesetup data onto the roms time vector
	int interpolate;
	int	time1, time2;
	double	y1;
	double	y2;
	double	mu;
	for(t=0;t<E->nTimeRoms;t++){
		// assume we are going to time interpolate
		interpolate = TRUE;
		time1 = -1;
		time2 = -1;

		// find the waveTimes that bound the current romsTime
		// get lower bound
		for(i=0;i<E->nTimeWaves;i++){
			if(E->waveTimes[i] <= E->romsTime[t])	// check this
				time1 = i;
		}
		// get upper bound
		for(i=E->nTimeWaves-1;i>=0;i--){
			if(E->waveTimes[i] >= E->romsTime[t])	// check this
				time2 = i;
		}
		// check to make sure that we need to interpolate
		// for some data, the time vectors are already aligned
		if(time1 == time2)
			interpolate = FALSE;

		if( (time1 < 0) || (time2 < 0)){
			printf("ERROR: wave time is not within the roms time vector!\n");
			printf("not extrapolating - exiting!\n");
			exit(1);
		}

		// setup the mu parameter for this time level
		// mu is the roms time value (between 0 and 1) on the waveTime interp_value
		//mu = (((OldValue - OldMin) * (NewMax - NewMin)) / (OldMax - OldMin)) + NewMin ;
		mu = ((E->romsTime[t] - E->waveTime[time1]) / (E->waveTime[time2] - E->waveTime[time1]));
		printf("lower bound (wave time) = %f\n", E->waveTime[time1] );
		printf("upper bound (wave time) = %f\n",  E->waveTime[time2]);
		printf("target roms time = %f\n", E->romsTime[t]);
		printf("mu = %f\n",mu);
		for(i=0;i<E->nLonRho;i++){
			for(j=0;j<E->nLatRho;j++){
				// interp the wave setup onto the roms time
				// y1 is the first vale at the wave time
				y1 = E->setup_on_roms[time1][i][j];
				// y2 is the second value at the wave time
				y1 = E->setup_on_roms[time2][i][j];


				if(interpolate == TRUE)
					//E->setup_on_roms_time_interp = LinearInterpolate( double y1,double y2, double mu);
				else
					E->setup_on_roms_time_interp = E->setup_on_roms;
			}
		}
	}




	// write out the derived coast line
	//E->nn_nx * E->nn_ny, E->nn_interp,
	int lat_dimid, lon_dimid, time_dimid, time_waves_dimid, dimIds[2], dimIds3d[3], dimIds3d_waves[3] ;
	int lat_varid, lon_varid, wave_time_varid, coastline_varid, setup_coast_varid, zeta_coast_varid;
	int ocean_time_varid;
	// create the file
	nc_create("roms_coast.nc", NC_CLOBBER, &ncid);
	// def dimensions
	nc_def_dim(ncid, "ocean_time", E->nTimeRoms, &time_dimid);
	nc_def_dim(ncid, "wave_time", E->nTimeWaves, &time_waves_dimid);
	nc_def_dim(ncid, "xi_rho", E->nLatRho, &lat_dimid);
	nc_def_dim(ncid, "eta_rho", E->nLonRho, &lon_dimid);
	// def vars
	dimIds[0] = lon_dimid;
	dimIds[1] = lat_dimid;

	dimIds3d[0] = time_dimid;
	dimIds3d[1] = lon_dimid;
	dimIds3d[2] = lat_dimid;

	dimIds3d_waves[0] = time_waves_dimid;
	dimIds3d_waves[1] = lon_dimid;
	dimIds3d_waves[2] = lat_dimid;

	nc_def_var(ncid, "wave_time", NC_DOUBLE, 1, &dimIds3d_waves[0], &wave_time_varid);
	nc_put_att_text(ncid, wave_time_varid, "long_name", strlen("time"), "time");
	nc_put_att_text(ncid, wave_time_varid, "units", strlen(E->roms_time_units), E->roms_time_units);
	nc_put_att_text(ncid, wave_time_varid, "calendar", strlen("gregorian"), "gregorian");

	nc_def_var(ncid, "ocean_time", NC_DOUBLE, 1, &dimIds3d[0], &ocean_time_varid);
	nc_put_att_text(ncid, ocean_time_varid, "long_name", strlen("time"), "time");
	nc_put_att_text(ncid, ocean_time_varid, "units", strlen(E->roms_time_units), E->roms_time_units);
	nc_put_att_text(ncid, ocean_time_varid, "calendar", strlen("gregorian"), "gregorian");

	nc_def_var(ncid, "lat_rho", NC_DOUBLE, 2, dimIds, &lat_varid);
	nc_put_att_text(ncid, lat_varid, "_CoordinateAxisType", strlen("Lat"), "Lat");
	nc_def_var(ncid, "lon_rho", NC_DOUBLE, 2, dimIds, &lon_varid);
	nc_put_att_text(ncid, lon_varid, "_CoordinateAxisType", strlen("Lon"), "Lon");
	nc_def_var(ncid, "coastline_mask", NC_DOUBLE, 2, dimIds, &coastline_varid);
	nc_put_att_text(ncid, coastline_varid, "coordinates", strlen("lat_rho lon_rho"), "lat_rho lon_rho");
	nc_def_var(ncid, "zeta_coast", NC_DOUBLE, 3, dimIds3d, &zeta_coast_varid);
	nc_put_att_text(ncid, zeta_coast_varid, "coordinates", strlen("lat_rho lon_rho"), "lat_rho lon_rho");
	nc_def_var(ncid, "setup", NC_DOUBLE, 3, dimIds3d_waves, &setup_coast_varid);
	nc_put_att_text(ncid, setup_coast_varid, "coordinates", strlen("lat_rho lon_rho"), "lat_rho lon_rho");
	nc_enddef(ncid);
	// write the data

	nc_put_var_double(ncid, ocean_time_varid, &E->romsTime[0]);
	nc_put_var_double(ncid, wave_time_varid, &E->wavesTime[0]);
	nc_put_var_double(ncid, lat_varid, &E->lat_rho[0][0]);
	nc_put_var_double(ncid, lon_varid, &E->lon_rho[0][0]);
	nc_put_var_double(ncid, coastline_varid, &E->coastline_mask[0][0]);
	nc_put_var_double(ncid, zeta_coast_varid, &zeta_coast[0][0][0]);
	nc_put_var_double(ncid, setup_coast_varid, &E->setup_on_roms[0][0][0]);

	// close the file
	nc_close(ncid);



	exit(1);




















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
	E->tide_time_units[attlen] = '\x0';
	printf("tide time units = %s\n", E->tide_time_units);

	// Make the Calendar calls
	// Parse the units strings
    if( (E->tide_ref_time = ut_parse( u_system, E->tide_time_units, UT_ASCII )) == NULL ) {
        fprintf( stderr, "Error parsing units string \"%s\"\n", E->roms_time_units );
        exit(-1);
        }

	// put the tide time on the same time units as the roms time
	for(t=0;t<E->nTimeTide;t++){

		//printf("PRE: tideTime[t] = %f  ", E->tideTime[t]);
		// convert tide time into year, month, day, hour, minute and seconds
	    if( (ierr = utCalendar2_cal( E->tideTime[t], E->tide_ref_time, &yr, &mo, &day, &hr, &min, &sec, calendar )) != 0 ) {
	        fprintf( stderr, "Error on utCalendar2_cal call: %d\n", ierr );
	        exit(-1);
	        }
	    //printf( "this date is %04d-%02d-%02d %02d:%02d:%06.3lf\n",yr, mo, day, hr, min, sec );

		// convert this date to be on the same units as the roms time
		if( (ierr = utInvCalendar2_cal( yr, mo, day, hr, min, sec, E->roms_ref_time, &E->tideTime[t], calendar )) != 0 ) {
			fprintf( stderr, "Error on utCalendar2_cal call: %d\n", ierr );
			exit(-1);
			}
		//printf( "POST: %04d-%02d-%02d %02d:%02d:%06.3lf is %lf %s in the %s calendar\n",
		//	yr, mo, day, hr, min, sec, E->tideTime[t], E->roms_time_units, calendar );
	}

	// find the index bounds where the tide time overlaps the roms time
	// the tide times should fully cover the roms times
	E->start_time_index = -1;
	E->start_time_roms = E->romsTime[0];

	// get the time start index for the tide file
	for(t=0;t<E->nTimeTide;t++){
		if(E->tideTime[t]<=E->start_time_roms)
			E->start_time_index = t;
	}
	if(E->start_time_index == -1){
		fprintf(stderr,"couldn't find a matching start time in the tide file.\n");
		fprintf(stderr,"check to make sure the tide file times sufficiently overlap\n");
		fprintf(stderr,"the ROMS times.\n\n");
		exit(1);
	}

	// get the end index for the tide file
	E->end_time_index = -1;
	E->end_time_roms = E->romsTime[E->nTimeRoms-1];
	for(t=E->nTimeTide-1;t>=0;t--){
		printf("t = %d, tide_time = %f\n",t,E->tideTime[t]);
		if(E->tideTime[t] >= E->end_time_roms)
			E->end_time_index = t;
	}

	if(E->end_time_index == -1){
		fprintf(stderr,"couldn't find a matching end time in the tide file.\n");
		fprintf(stderr,"check to make sure the tide file times sufficiently overlap\n");
		fprintf(stderr,"the ROMS times.\n\n");
		fprintf(stderr,"end time ROMS = %f\n", E->end_time_roms);
		fprintf(stderr,"end time tide = %f\n", E->tideTime[E->nTimeTide-1]);
		exit(1);
	}

	printf("start index = %d\n", E->start_time_index);
	printf("end index = %d\n", E->end_time_index);

	// make the time vector for the output file
	E->interp_time = malloc((E->end_time_index-E->start_time_index)*sizeof(double));
	count=0;
	for(t=E->start_time_index;t<E->end_time_index;t++){
		E->interp_time[count] = E->tideTime[t];
		count++;
	}


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

	// malloc memory for the lib-nn arrays and the loop working arrays
	E->nn_diff = malloc(E->nLonTide * E->nLatTide * sizeof(point));
	// and malloc room for the interpolated data
	// we already know the mesh dimensions from the input tide file
	// E->nx = number of longitudes in file
	// E->ny = number of latitudes in file
	E->nn_interp = malloc(E->nLonTide * E->nLatTide * sizeof(point));

	E->tide_on_roms = malloc3d_double((E->end_time_index-E->start_time_index), E->nLonRho, E->nLatRho);

	printf("(E->end_time_index-E->start_time_index) = %d\n", E->end_time_index-E->start_time_index);


	// set up variables for the lib-nn calls
	// now do the interpolation via lib-nn
	E->nn_weight = 0.0;
	// so have E->weight = 0.0; somewhere
	// E->diff is the input data list
	// E->interp is the output data list
	// E-> nx and E->ny is the grid dimensions
	E->nn_dx = fabs(E->tideLat[1] - E->tideLat[0]); // lon
	E->nn_dy = fabs(E->tideLon[1] - E->tideLon[0]); // lat

	// for each time level
	for(t=0;t<(E->end_time_index-E->start_time_index);t++){
		// find out how many valid data points we have
		// and setup the input array for lib-nn
		E->nn_n = 0;
		for(i=0;i<E->nLatTide;i++){
			for(j=0;j<E->nLonTide;j++){
				if(E->tide_data[t][0][i][j] > -999.0){
					E->nn_diff[E->nn_n].y = E->tideLat[i];
					E->nn_diff[E->nn_n].x = E->tideLon[j];
					E->nn_diff[E->nn_n].z = E->tide_data[t+E->start_time_index][0][i][j];

					E->nn_n++;
				}
			}
		}

		// do the natural neighbour interpolation
		// this will interpolate our scattered point data
		// onto an orthogonal grid with dimensions
		// E->nx x E->ny
		// figure out E->nx and E->ny
		get_mesh_dimensions(E, E->nn_n, E->nn_diff);
		// nn_interp is the interpolated data from lib-nn
		nn_interp = malloc(E->nn_nx * E->nn_ny * sizeof(double));
		printf("nn interp field: nn_nx = %d, nn_ny = %d\n", E->nn_nx, E->nn_ny);


		printf("doing nn interpolation...\n"); fflush(stdout);
		// E->n is the number of points in E->diff
		// E->weight is a constant (equal to 0.0)

		nn_interp_to_mesh(E, E->nn_n, E->nn_weight, E->nn_diff, E->nLonTide, E->nLatTide, E->nn_interp);
		printf("done\n"); fflush(stdout);


		// put the interpolated data into the tide array so we can write it out to
		// a netcdf file
		for(i=0;i<E->nn_nx*E->nn_ny;i++){
			//fprintf(fptr,"%.15g %.15g %.15g\n", p[i].x, p[i].y, p[i].z);
			nn_interp[i] = E->nn_interp[i].z;
		}

		// interpolate the tide data for each lon_rho and lat_rho point
		printf("t = %d\n", t);

		for(i=0;i<E->nLonRho;i++)
			for(j=0;j<E->nLatRho;j++)
				E->tide_on_roms[t][i][j] = 0.0;


		// temporarily splat the nn_interped field over the original data
		count = 0;
		for(i=0;i<E->nLatTide;i++){
			for(j=0;j<E->nLonTide;j++){
				E->tide_data[t][0][i][j] = nn_interp[count];
				count++;
			}
		}

		// interpolate tide to roms grid
	    interp_tide_to_roms(E,t);


		// apply land-sea mask
		for(i=0;i<E->nLonRho;i++){
			for(j=0;j<E->nLatRho;j++){
				if(E->mask_rho[i][j] == 0)
					E->tide_on_roms[t][i][j] = NC_FILL_DOUBLE;
			}
		}
		free(nn_interp);

	} // end of loop over tide time levels

	// write the interpolated field to file
	write_netcdf(E);


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
