#include "grid.h"

void process_tides(e *E){

    int	ncid;
	int varid;
	int retval;
	size_t attlen = 0;
    size_t from[4];
    size_t to[4];

    double	sec;
    int	ierr, yr, mo, day, hr, min;

    // nn variables
	double *nn_interp;

    int     i,j,t;
	int count;

    static int  beenHere = FALSE;


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

    // process the spatial dimensions
    E->tideLon = malloc(E->nLonTide*sizeof(double));
	E->tideLat = malloc(E->nLatTide*sizeof(double));

    nc_inq_varid(ncid, "lat", &varid);
    if((retval = nc_get_var_double(ncid, varid, &E->tideLat[0])))
		fail("failed to read tide lat data data: error is %d\n",retval);

	nc_inq_varid(ncid, "lon", &varid);
    if((retval = nc_get_var_double(ncid, varid, &E->tideLon[0])))
		fail("failed to read tide lon data data: error is %d\n",retval);

    // get spatial bounds
    // find the dimension indexes that cover the spatial region we need
    int lat_start;
    for(i=0;i<E->nLatTide;i++){
        if(E->tideLat[i] < (E->roms_min_lat-0.15))
            lat_start = i;
    }
    int lat_end;
    for(i=E->nLatTide-1;i>=0;i--){
        if(E->tideLat[i] > (E->roms_max_lat+0.15))
            lat_end = i;
    }
    printf("tide data start lat = %f (%d), end lat = %f (%d)\n", E->tideLat[lat_start],lat_start, E->tideLat[lat_end],lat_end);
    int lon_start;
    for(i=0;i<E->nLonTide;i++){
        if(E->tideLon[i] < (E->roms_min_lon-0.15))
            lon_start = i;
    }
    int lon_end;
    for(i=E->nLonTide-1;i>=0;i--){
        if(E->tideLon[i] > (E->roms_max_lon+0.15))
            lon_end = i;
    }

    printf("tide data start lon = %f, end lon = %f\n", E->tideLon[lon_start], E->tideLon[lon_end]);

    // TODO: add some error checking to the bounds code.
    // for example, if the spatial extent does not overlap then throw an exception

    // now just read in what we want from the files
    free(E->tideLat);
    free(E->tideLon);
    E->nLatTide = (lat_end - lat_start);
    E->nLonTide = (lon_end - lon_start);
    E->tideLat = malloc(E->nLatTide*sizeof(double));
    E->tideLon = malloc(E->nLonTide*sizeof(double));
    // now re-read the lat and lon data
    size_t spatial_from[1], spatial_to[1];
    spatial_from[0] = lat_start;    spatial_to[0] = lat_end-lat_start;
    nc_inq_varid(ncid, "lat", &varid);
    if((retval = nc_get_vara_double(ncid, varid, spatial_from, spatial_to, &E->tideLat[0])))
        fail("failed to read tides lat data: error is %d\n", retval);
    spatial_from[0] = lon_start;    spatial_to[0] = lon_end-lon_start;
    nc_inq_varid(ncid, "lon", &varid);
    if((retval = nc_get_vara_double(ncid, varid, spatial_from, spatial_to, &E->tideLon[0])))
        fail("failed to read tides lon data: error is %d\n", retval);



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
    if( (E->tide_ref_time = ut_parse( E->u_system, E->tide_time_units, UT_ASCII )) == NULL ) {
        fprintf( stderr, "Error parsing units string \"%s\"\n", E->roms_time_units );
        exit(-1);
        }

	// put the tide time on the same time units as the roms time
	for(t=0;t<E->nTimeTide;t++){

		//printf("PRE: tideTime[t] = %f  ", E->tideTime[t]);
		// convert tide time into year, month, day, hour, minute and seconds
	    if( (ierr = utCalendar2_cal( E->tideTime[t], E->tide_ref_time, &yr, &mo, &day, &hr, &min, &sec, E->calendar )) != 0 ) {
	        fprintf( stderr, "Error on utCalendar2_cal call: %d\n", ierr );
	        exit(-1);
	        }
	    //printf( "this date is %04d-%02d-%02d %02d:%02d:%06.3lf\n",yr, mo, day, hr, min, sec );

		// convert this date to be on the same units as the roms time
		if( (ierr = utInvCalendar2_cal( yr, mo, day, hr, min, sec, E->roms_ref_time, &E->tideTime[t], E->calendar )) != 0 ) {
			fprintf( stderr, "Error on utCalendar2_cal call: %d\n", ierr );
			exit(-1);
			}
		//printf( "POST: %04d-%02d-%02d %02d:%02d:%06.3lf is %lf %s in the %s calendar\n",
		//	yr, mo, day, hr, min, sec, E->tideTime[t], E->roms_time_units, calendar );
	}

	// find the index bounds where the tide time overlaps the roms time
	// the tide times should fully cover the roms times
	E->tide_start_time_index = -1;
	E->start_time_roms = E->romsTime[0];

	// get the time start index for the tide file
	for(t=0;t<E->nTimeTide;t++){
		if(E->tideTime[t]<=E->start_time_roms)
			E->tide_start_time_index = t;
	}
	if(E->tide_start_time_index == -1){
		fprintf(stderr,"couldn't find a matching start time in the tide file.\n");
		fprintf(stderr,"check to make sure the tide file times sufficiently overlap\n");
		fprintf(stderr,"the ROMS times.\n\n");
		exit(1);
	}

	// get the end index for the tide file
	E->tide_end_time_index = -1;
	E->end_time_roms = E->romsTime[E->nTimeRoms-1];
	for(t=E->nTimeTide-1;t>=0;t--){
		//printf("t = %d, tide_time = %f\n",t,E->tideTime[t]);
		if(E->tideTime[t] >= E->end_time_roms)
			E->tide_end_time_index = t;
	}

	if(E->tide_end_time_index == -1){
		fprintf(stderr,"couldn't find a matching end time in the tide file.\n");
		fprintf(stderr,"check to make sure the tide file times sufficiently overlap\n");
		fprintf(stderr,"the ROMS times.\n\n");
		fprintf(stderr,"end time ROMS = %f\n", E->end_time_roms);
		fprintf(stderr,"end time tide = %f\n", E->tideTime[E->nTimeTide-1]);
		exit(1);
	}

	printf("Tide start index = %d\n", E->tide_start_time_index);
	printf("Tide end index = %d\n", E->tide_end_time_index);

    E->nTimeTideSubset = (E->tide_end_time_index-E->tide_start_time_index+1);

	// make the time vector for the output file
	E->tide_interp_time = malloc(E->nTimeTideSubset*sizeof(double));
	count=0;
	for(t=E->tide_start_time_index;t<E->tide_end_time_index;t++){
		E->tide_interp_time[count] = E->tideTime[t];
		count++;
	}


    // read the data
    E->tide_data = malloc4d_double(E->nTimeTideSubset, E->nLevTide, E->nLatTide, E->nLonTide);

    from[0] = E->tide_start_time_index;     to[0] = E->nTimeTideSubset;
    from[1] = 0;                            to[1] = E->nLevTide;
    from[2] = 0;                            to[2] = E->nLatTide;
    from[3] = 0;                            to[3] = E->nLonTide;

    nc_inq_varid(ncid, "z", &varid);
    // should only read in the data we need
    if((retval = nc_get_vara_double(ncid, varid, from, to, &E->tide_data[0][0][0][0])))
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

	E->tide_on_roms = malloc3d_double(E->nTimeTideSubset, E->nLonRho, E->nLatRho);

	printf("(E->end_time_index-E->start_time_index) = %d\n", E->nTimeTideSubset);


	// set up variables for the lib-nn calls
	// now do the interpolation via lib-nn
	E->nn_weight = 0.0;
	// so have E->weight = 0.0; somewhere
	// E->diff is the input data list
	// E->interp is the output data list
	// E-> nx and E->ny is the grid dimensions
	E->nn_dx = fabs(E->tideLat[1] - E->tideLat[0]); // lon
	E->nn_dy = fabs(E->tideLon[1] - E->tideLon[0]); // lat

    E->nn_rule = NON_SIBSONIAN; //SIBSON;

	// for each time level
	for(t=0;t<E->nTimeTideSubset;t++){
		// find out how many valid data points we have
		// and setup the input array for lib-nn
		E->nn_n = 0;
		for(i=0;i<E->nLatTide;i++){
			for(j=0;j<E->nLonTide;j++){
				if(E->tide_data[t][0][i][j] > -999.0){
					E->nn_diff[E->nn_n].y = E->tideLat[i];
					E->nn_diff[E->nn_n].x = E->tideLon[j];
					E->nn_diff[E->nn_n].z = E->tide_data[t][0][i][j];

					E->nn_n++;
				}
			}
		}

		// do the natural neighbour interpolation
		// this will interpolate our scattered point data
		// onto an orthogonal grid with dimensions
		// E->nx x E->ny

		// E->nn_interp is the interpolated data from lib-nn
        if(beenHere == FALSE){
            // figure out E->nx and E->ny
    		get_mesh_dimensions(E, E->nn_n, E->nn_diff);
		    printf("nn interp field: nn_nx = %d, nn_ny = %d\n", E->nn_nx, E->nn_ny);
            beenHere = TRUE;
        }

		printf("doing nn interpolation...\n"); fflush(stdout);
		// E->n is the number of points in E->diff
		// E->weight is a constant (equal to 0.0)

		nn_interp_to_mesh(E, E->nn_n, E->nn_weight, E->nn_diff, E->nLonTide, E->nLatTide, E->nn_interp);
		printf("done\n"); fflush(stdout);


		// put the interpolated data into the tide array so we can write it out to
		// a netcdf file
		// interpolate the tide data for each lon_rho and lat_rho point
		printf("t = %d\n", t);

		// temporarily splat the nn_interped field over the original data
		count = 0;
		for(i=0;i<E->nLatTide;i++){
			for(j=0;j<E->nLonTide;j++){
				E->tide_data[t][0][i][j] = E->nn_interp[count].z;
				count++;
			}
		}

		// interpolate tide to roms grid
	    interp_tide_to_roms(E,t);

        // apply coastline mask to tide data
        for(i=0;i<E->nLonRho;i++){
			for(j=0;j<E->nLatRho;j++){
				if(E->coastline_mask[i][j]  == 0)
					E->tide_on_roms[t][i][j] = NC_FILL_DOUBLE;
			}
		}



	} // end of loop over tide time levels
    // time interp tide data to roms time

    free(E->tideLon);
    free(E->tideLat);
    free(E->tide_data);
    free(E->nn_interp);
    free(E->nn_diff);
}
