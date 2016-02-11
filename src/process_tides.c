#include "grid.h"

void process_tides(e *E){

    int	ncid;
	int varid;
	int retval;
	size_t attlen = 0;

    double	sec;
    int	ierr, yr, mo, day, hr, min;

    // nn variables
	double *nn_interp;

    int     i,j,t;
	int count;


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
		printf("t = %d, tide_time = %f\n",t,E->tideTime[t]);
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

	printf("start index = %d\n", E->tide_start_time_index);
	printf("end index = %d\n", E->tide_end_time_index);

	// make the time vector for the output file
	E->tide_interp_time = malloc((E->tide_end_time_index-E->tide_start_time_index)*sizeof(double));
	count=0;
	for(t=E->tide_start_time_index;t<E->tide_end_time_index;t++){
		E->tide_interp_time[count] = E->tideTime[t];
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

	E->tide_on_roms = malloc3d_double((E->tide_end_time_index-E->tide_start_time_index), E->nLonRho, E->nLatRho);

	printf("(E->end_time_index-E->start_time_index) = %d\n", E->tide_end_time_index-E->tide_start_time_index);


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
	for(t=0;t<(E->tide_end_time_index-E->tide_start_time_index);t++){
		// find out how many valid data points we have
		// and setup the input array for lib-nn
		E->nn_n = 0;
		for(i=0;i<E->nLatTide;i++){
			for(j=0;j<E->nLonTide;j++){
				if(E->tide_data[t][0][i][j] > -999.0){
					E->nn_diff[E->nn_n].y = E->tideLat[i];
					E->nn_diff[E->nn_n].x = E->tideLon[j];
					E->nn_diff[E->nn_n].z = E->tide_data[t+E->tide_start_time_index][0][i][j];

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

}
