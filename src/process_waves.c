#include "grid.h"

void process_auswave(e *E){

    // netcdf vars
    int	ncid;
	int varid;
	int retval;
	size_t attlen = 0;
    size_t from[3];
    size_t to[3];

    // for time conversion
	static	char *calendar = "Standard";
    ut_system	*u_system;
    double	sec;
    int	ierr, yr, mo, day, hr, min;

    // nn variables
	double *nn_interp;

    int     i,j,t;
	int count;

    static int  beenHere = FALSE;


    // read in the auswave data so we can estimate wave setup
	//time = UNLIMITED ; // (192 currently)
	//lat = 411 ;
	//lon = 441 ;
	// open the file
    if((retval = nc_open(E->wave_input, NC_NOWRITE, &ncid)))
        fail("failed to open wave setup input file: error is %d\n",retval);

	// get the time data
	if((retval = nc_inq_dimid(ncid, "time", &varid)))
        fail("failed to get auswave dimid: error is %d\n",retval);

    if((retval = nc_inq_dimlen(ncid,varid,&E->nTimeWaves)))
        fail("failed to get auswave lat dimlen: error is %d\n",retval);

    printf("aus waves_times = %zu\n", E->nTimeWaves);

	// get the auswave lat data
	if((retval = nc_inq_dimid(ncid, "lat", &varid)))
        fail("failed to get auswave lat dimid: error is %d\n",retval);

    if((retval = nc_inq_dimlen(ncid,varid,&E->nLatWaves)))
        fail("failed to get auswave lat dimlen: error is %d\n",retval);

    printf("auswave lat = %zu\n", E->nLatWaves);

	// get the auswave lon data
	if((retval = nc_inq_dimid(ncid, "lon", &varid)))
        fail("failed to get auswave lon dimid: error is %d\n",retval);

    if((retval = nc_inq_dimlen(ncid,varid,&E->nLonWaves)))
        fail("failed to get auswave lon dimlen: error is %d\n",retval);

    printf("auswave lat = %zu\n", E->nLonWaves);

    // process spatial dimensions
    // malloc room for the dim variable arrays
	E->wavesLat = malloc(E->nLatWaves*sizeof(double));
	E->wavesLon = malloc(E->nLonWaves*sizeof(double));

	nc_inq_varid(ncid, "lat", &varid);
    if((retval = nc_get_var_double(ncid, varid, &E->wavesLat[0])))
		fail("failed to read waves lat data: error is %d\n", retval);

		//printf("waves lat[0] = %f\n", E->wavesLat[0]);
    nc_inq_varid(ncid, "lon", &varid);
    if((retval = nc_get_var_double(ncid, varid, &E->wavesLon[0])))
		fail("failed to read waves lon data: error is %d\n", retval);

	//printf("waves lon[0] = %f\n", E->wavesLon[0]);

    // find the dimension indexes that cover the spatial region we need
    int lat_end;
    for(i=0;i<E->nLatWaves;i++){
        if(E->wavesLat[i] > (E->roms_min_lat-0.15))    // greater than because auswave lat is monotonically decreasing
            lat_end = i;
    }
    int lat_start;
    for(i=E->nLatWaves-1;i>=0;i--){
        if(E->wavesLat[i] < (E->roms_max_lat+0.15))    // less than because auswave lat is monotonically decreasing
            lat_start = i;
    }
    printf("wave data start lat = %f (%d), end lat = %f (%d)\n", E->wavesLat[lat_start],lat_start, E->wavesLat[lat_end],lat_end);
    int lon_start;
    for(i=0;i<E->nLonWaves;i++){
        if(E->wavesLon[i] < (E->roms_min_lon-0.15))
            lon_start = i;
    }
    int lon_end;
    for(i=E->nLonWaves-1;i>=0;i--){
        if(E->wavesLon[i] > (E->roms_max_lon+0.15))
            lon_end = i;
    }

    printf("wave data start lon = %f, end lon = %f\n", E->wavesLon[lon_start], E->wavesLon[lon_end]);

    // TODO: add some error checking to the bounds code.
    // for example, if the spatial extent does not overlap then throw an exception

    // now just read in what we want from the files
    free(E->wavesLat);
    free(E->wavesLon);
    E->nLatWaves = (lat_end - lat_start);
    E->nLonWaves = (lon_end - lon_start);
    E->wavesLat = malloc(E->nLatWaves*sizeof(double));
    E->wavesLon = malloc(E->nLonWaves*sizeof(double));
    // now re-read the lat and lon data
    size_t spatial_from[1], spatial_to[1];
    spatial_from[0] = lat_start;    spatial_to[0] = lat_end-lat_start;
    nc_inq_varid(ncid, "lat", &varid);
    if((retval = nc_get_vara_double(ncid, varid, spatial_from, spatial_to, &E->wavesLat[0])))
        fail("failed to read waves lat data: error is %d\n", retval);
    spatial_from[0] = lon_start;    spatial_to[0] = lon_end-lon_start;
    nc_inq_varid(ncid, "lon", &varid);
    if((retval = nc_get_vara_double(ncid, varid, spatial_from, spatial_to, &E->wavesLon[0])))
        fail("failed to read waves lon data: error is %d\n", retval);

    // process time
	E->wavesTime = malloc(E->nTimeWaves*sizeof(double));
    // read the data from the waves output file
	nc_inq_varid(ncid, "time", &varid);
    if((retval = nc_get_var_double(ncid, varid, &E->wavesTime[0])))
		fail("failed to read auswave time data: error is %d\n", retval);

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
    if( (E->waves_ref_time = ut_parse( E->u_system, E->waves_time_units, UT_ASCII )) == NULL ) {
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

	// put the waves time on the same time units as the roms time
	for(t=0;t<E->nTimeWaves;t++){

		//printf("PRE: waveTimes[t] = %f  ", E->wavesTime[t]);
		// convert tide time into year, month, day, hour, minute and seconds
	    if( (ierr = utCalendar2_cal( E->wavesTime[t], E->waves_ref_time, &yr, &mo, &day, &hr, &min, &sec, E->calendar )) != 0 ) {
	        fprintf( stderr, "Error on utCalendar2_cal call: %d\n", ierr );
	        exit(-1);
	        }
	    //printf( "this date is %04d-%02d-%02d %02d:%02d:%06.3lf\n",yr, mo, day, hr, min, sec );

		// convert this date to be on the same units as the roms time
		if( (ierr = utInvCalendar2_cal( yr, mo, day, hr, min, sec, E->roms_ref_time, &E->wavesTime[t], E->calendar )) != 0 ) {
			fprintf( stderr, "Error on utCalendar2_cal call: %d\n", ierr );
			exit(-1);
			}
		//printf( "POST: %04d-%02d-%02d %02d:%02d:%06.3lf is %lf %s in the %s calendar\n",
		//	yr, mo, day, hr, min, sec, E->wavesTime[t], E->roms_time_units, E->calendar );
	}


    // find the index bounds where the waves time overlaps the roms time
	// the waves times should fully cover the roms times
	E->waves_start_time_index = -1;
	E->start_time_roms = E->romsTime[0];

	// get the time start index for the tide file
	for(t=0;t<E->nTimeWaves;t++){
		if(E->wavesTime[t]<=E->start_time_roms)
			E->waves_start_time_index = t;
	}
	if(E->waves_start_time_index == -1){
		fprintf(stderr,"couldn't find a matching start time in the waves file.\n");
		fprintf(stderr,"check to make sure the waves file times sufficiently overlap\n");
		fprintf(stderr,"the ROMS times.\n\n");
		exit(1);
	}

	// get the end index for the tide file
	E->waves_end_time_index = -1;
	E->end_time_roms = E->romsTime[E->nTimeRoms-1];
	for(t=E->nTimeWaves-1;t>=0;t--){
		//printf("t = %d, wave_time = %f\n",t,E->wavesTime[t]);
		if(E->wavesTime[t] >= E->end_time_roms)
			E->waves_end_time_index = t;
	}

	if(E->waves_end_time_index == -1){
		fprintf(stderr,"couldn't find a matching end time in the waves file.\n");
		fprintf(stderr,"check to make sure the wave file times sufficiently overlap\n");
		fprintf(stderr,"the ROMS times.\n\n");
		fprintf(stderr,"end time ROMS = %f\n", E->end_time_roms);
		fprintf(stderr,"end time waves = %f\n", E->wavesTime[E->nTimeWaves-1]);
		exit(1);
	}

	printf("start index = %d\n", E->waves_start_time_index);
	printf("end index = %d\n", E->waves_end_time_index);
    E->nTimeWavesSubset = (E->waves_end_time_index - E->waves_start_time_index)+1;
    // malloc enough room for the variable arrays
    //sig_wav_ht(time, lat, lon)
	E->Hs = malloc3d_double(E->nTimeWavesSubset, E->nLatWaves, E->nLonWaves);
	E->Tp = malloc3d_double(E->nTimeWavesSubset, E->nLatWaves, E->nLonWaves);

	// make the time vector for the output file
	E->waves_interp_time = malloc(E->nTimeWavesSubset*sizeof(double));
	count=0;
	for(t=E->waves_start_time_index;t<=E->waves_end_time_index;t++){
		E->waves_interp_time[count] = E->wavesTime[t];
		count++;
	}


	// get the sig wave height
	nc_inq_varid(ncid, "sig_wav_ht", &varid);
    //if((retval = nc_get_var_double(ncid, varid, &E->Hs[0][0][0])))
	//	fail("failed to read waves setup data: error is %d\n", retval);

    from[0] = E->waves_start_time_index;    to[0] = E->nTimeWavesSubset;
    from[1] = 0;                            to[1] = E->nLatWaves;
    from[2] = 0;                            to[2] = E->nLonWaves;
    if((retval = nc_get_vara_double(ncid, varid, from, to, &E->Hs[0][0][0])))
        fail("failed to read waves Hs data: error is %d\n", retval);

	printf("sig_wave_ht[0][0][0] = %f\n", E->Hs[0][0][0]);

	// get the peak period
	nc_inq_varid(ncid, "pk_wav_per", &varid);
    if((retval = nc_get_vara_double(ncid, varid, from, to, &E->Tp[0][0][0])))
        fail("failed to read waves Hs data: error is %d\n", retval);

	printf("pk_wav_per[0][0][0] = %f\n", E->Tp[0][0][0]);

	// close the file
	nc_close(ncid);

    // flip the auswave data so the lat vector is monotonically increasing
    double ***flipData = malloc3d_double(E->nTimeWavesSubset, E->nLatWaves, E->nLonWaves);
    double  *flipLat = malloc(E->nLatWaves*sizeof(double));
    // flip the lat vector
    for(i=0;i<E->nLatWaves;i++){
        flipLat[i] = E->wavesLat[E->nLatWaves-1-i];
    }
    // copy the flipped data back
    for(i=0;i<E->nLatWaves;i++){
        E->wavesLat[i] = flipLat[i];
    }
    // flip the Hs data array
    for(t=0;t<E->nTimeWavesSubset;t++){
        for(i=0;i<E->nLatWaves;i++){
            for(j=0;j<E->nLonWaves;j++){
                flipData[t][i][j] = E->Hs[t][E->nLatWaves-1-i][j];
            }
        }
    }
    // copy it back
    for(t=0;t<E->nTimeWavesSubset;t++){
        for(i=0;i<E->nLatWaves;i++){
            for(j=0;j<E->nLonWaves;j++){
                E->Hs[t][i][j] = flipData[t][i][j];
            }
        }
    }
    // flip the Tp data array
    for(t=0;t<E->nTimeWavesSubset;t++){
        for(i=0;i<E->nLatWaves;i++){
            for(j=0;j<E->nLonWaves;j++){
                flipData[t][i][j] = E->Tp[t][E->nLatWaves-1-i][j];
            }
        }
    }
    // copy it back
    for(t=0;t<E->nTimeWavesSubset;t++){
        for(i=0;i<E->nLatWaves;i++){
            for(j=0;j<E->nLonWaves;j++){
                E->Tp[t][i][j] = flipData[t][i][j];
            }
        }
    }
    free(flipData);
    free(flipLat);

    // malloc room for the output nearest neighbor interp auswave  data
	// target grid for auswave data data

    // do a natural neighbour interpolation on the tide data we just read in
    // to fill in the land masked values before interpolating onto the ROMS grid

    // malloc memory for the lib-nn arrays and the loop working arrays
    E->nn_diff = malloc(E->nLonWaves * E->nLatWaves * sizeof(point));
    // and malloc room for the interpolated data
    // we already know the mesh dimensions from the input tide file
    // E->nx = number of longitudes in file
    // E->ny = number of latitudes in file
    E->nn_interp = malloc(E->nLonWaves * E->nLatWaves * sizeof(point));

    E->Hs_on_roms = malloc3d_double(E->nTimeWavesSubset, E->nLonRho, E->nLatRho);
    E->Tp_on_roms = malloc3d_double(E->nTimeWavesSubset, E->nLonRho, E->nLatRho);
    // just for writing the netcdf file - trash this!
    //double *time_vector = malloc(E->nTimeWavesSubset*sizeof(double));
    printf("(E->waves_end_time_index-E->waves_start_time_index) = %d\n", E->nTimeWavesSubset);


    // set up variables for the lib-nn calls
    // now do the Hs interpolation via lib-nn
    E->nn_weight = 0.0;
    // so have E->weight = 0.0; somewhere
    // E->diff is the input data list
    // E->interp is the output data list
    // E-> nx and E->ny is the grid dimensions
    E->nn_dx = fabs(E->wavesLat[1] - E->wavesLat[0]); // lon
    E->nn_dy = fabs(E->wavesLon[1] - E->wavesLon[0]); // lat

    E->nn_rule = NON_SIBSONIAN;
    //printf("original data at i=74, j=237 Hs = %f\n", E->Hs[0][74][237]);
    //printf("nn_dx = %f, nn_dy = %f\n", E->nn_dx, E->nn_dy);

    // for each time level
    for(t=0;t<E->nTimeWavesSubset;t++){
        // find out how many valid data points we have
        // and setup the input array for lib-nn
        //time_vector[t] = (double)t;
        E->nn_n = 0;
        for(i=0;i<E->nLatWaves;i++){
            for(j=0;j<E->nLonWaves;j++){
                if(E->Hs[t][i][j] > -999.0){
                    E->nn_diff[E->nn_n].y = E->wavesLat[i];
                    E->nn_diff[E->nn_n].x = E->wavesLon[j];
                    E->nn_diff[E->nn_n].z = E->Hs[t][i][j];
                    //printf("i = %d, j = %d, lat = %.15g lon = %.15g Hs = %.15g\n", E->nn_diff[E->nn_n].x, E->nn_diff[E->nn_n].y, E->nn_diff[E->nn_n].z);
                    E->nn_n++;
                }
            }
        }

        // do the natural neighbour interpolation
        // this will interpolate our scattered point data
        // onto an orthogonal grid with dimensions
        // E->nx x E->ny
        // figure out E->nx and E->ny

        if(beenHere == FALSE){
            get_mesh_dimensions(E, E->nn_n, E->nn_diff);
            // nn_interp is the interpolated data from lib-nn
            beenHere = TRUE;
        }

        printf("doing nn interpolation (waves Hs)...\n"); fflush(stdout);
        // E->n is the number of points in E->diff
        // E->weight is a constant (equal to 0.0)

        nn_interp_to_mesh(E, E->nn_n, E->nn_weight, E->nn_diff, E->nLonWaves, E->nLatWaves, E->nn_interp);
        printf("done\n"); fflush(stdout);

        // interpolate the tide data for each lon_rho and lat_rho point
        printf("t = %d\n", t);

        // temporarily splat the nn_interped field over the original data
        count = 0;
        for(i=0;i<E->nLatWaves;i++){
            for(j=0;j<E->nLonWaves;j++){
                E->Hs[t][i][j] = E->nn_interp[count].z;
                count++;
            }
        }

        // interpolate Hs to roms grid
        interp_hs_to_roms(E,t);

        // apply land-sea mask
        for(i=0;i<E->nLonRho;i++){
            for(j=0;j<E->nLatRho;j++){
                if(E->mask_rho[i][j] == 0)
                    E->Hs_on_roms[t][i][j] = NC_FILL_DOUBLE;
            }
        }

        /*
        // write it out to check
        //E->nn_nx * E->nn_ny, E->nn_interp,
        int lat_dimid, lon_dimid, time_dimid, dimIds[2];
        int lat_varid, lon_varid, time_varid, interp_varid;
        // create the file
        nc_create("hs_interp.nc", NC_CLOBBER, &ncid);
        // def dimensions
        nc_def_dim(ncid, "lat", E->nn_ny-1, &lat_dimid);
        nc_def_dim(ncid, "lon", E->nn_nx-1, &lon_dimid);
        // def vars
        dimIds[0] = lat_dimid;
        dimIds[1] = lon_dimid;
        nc_def_var(ncid, "lat", NC_DOUBLE, 1, &dimIds[0], &lat_varid);
        nc_def_var(ncid, "lon", NC_DOUBLE, 1, &dimIds[1], &lon_varid);
        nc_def_var(ncid, "hs_interp", NC_DOUBLE, 2, dimIds, &interp_varid);
        nc_enddef(ncid);
        // write the data
        nc_put_var_double(ncid, lat_varid, &E->wavesLat[0]);
        nc_put_var_double(ncid, lon_varid, &E->wavesLon[0]);
        nc_put_var_double(ncid, interp_varid, nn_interp);
        // close the file
        nc_close(ncid);
        //exit(1);
        */
        //free(nn_interp);

    } // end of loop over Hs time levels


    // now interp the Tp variable
    // set up variables for the lib-nn calls
    // now do the interpolation via lib-nn
    E->nn_weight = 0.0;
    // so have E->weight = 0.0; somewhere
    // E->diff is the input data list
    // E->interp is the output data list
    // E-> nx and E->ny is the grid dimensions
    E->nn_dx = fabs(E->wavesLat[1] - E->wavesLat[0]); // lon
    E->nn_dy = fabs(E->wavesLon[1] - E->wavesLon[0]); // lat

    E->nn_rule = NON_SIBSONIAN;

    // for each time level
    for(t=0;t<E->nTimeWavesSubset;t++){
        // find out how many valid data points we have
        // and setup the input array for lib-nn
        //time_vector[t] = (double)t;
        E->nn_n = 0;
        for(i=0;i<E->nLatWaves;i++){
            for(j=0;j<E->nLonWaves;j++){
                if(E->Tp[t][i][j] > -999.0){
                    E->nn_diff[E->nn_n].y = E->wavesLat[i];
                    E->nn_diff[E->nn_n].x = E->wavesLon[j];
                    E->nn_diff[E->nn_n].z = E->Tp[t][i][j];

                    E->nn_n++;
                }
            }
        }

        // do the natural neighbour interpolation
        // this will interpolate our scattered point data
        // onto an orthogonal grid with dimensions
        // E->nx x E->ny
        // figure out E->nx and E->ny
        if(beenHere == FALSE){
            get_mesh_dimensions(E, E->nn_n, E->nn_diff);
            printf("nn interp field (waves Tp): nn_nx = %d, nn_ny = %d\n", E->nn_nx, E->nn_ny);
        }

        printf("doing nn interpolation (waves Tp)...\n"); fflush(stdout);
        // E->n is the number of points in E->diff
        // E->weight is a constant (equal to 0.0)

        nn_interp_to_mesh(E, E->nn_n, E->nn_weight, E->nn_diff, E->nLonWaves, E->nLatWaves, E->nn_interp);
        printf("done\n"); fflush(stdout);

        // interpolate the tide data for each lon_rho and lat_rho point
        printf("t = %d\n", t);

        // temporarily splat the nn_interped field over the original data
        count = 0;
        for(i=0;i<E->nLatWaves;i++){
            for(j=0;j<E->nLonWaves;j++){
                E->Tp[t][i][j] = E->nn_interp[count].z;
                count++;
            }
        }

        // interpolate tp to roms grid
        interp_tp_to_roms(E,t);

        // apply land-sea mask
        for(i=0;i<E->nLonRho;i++){
            for(j=0;j<E->nLatRho;j++){
                if(E->mask_rho[i][j] == 0)
                    E->Tp_on_roms[t][i][j] = NC_FILL_DOUBLE;
            }
        }
        /*
        // write it out to check
        //E->nn_nx * E->nn_ny, E->nn_interp,
        int lat_dimid, lon_dimid, time_dimid, dimIds[2];
        int lat_varid, lon_varid, time_varid, interp_varid;
        // create the file
        nc_create("tp_interp.nc", NC_CLOBBER, &ncid);
        // def dimensions
        nc_def_dim(ncid, "lat", E->nn_ny-1, &lat_dimid);
        nc_def_dim(ncid, "lon", E->nn_nx-1, &lon_dimid);
        // def vars
        dimIds[0] = lat_dimid;
        dimIds[1] = lon_dimid;
        nc_def_var(ncid, "lat", NC_DOUBLE, 1, &dimIds[0], &lat_varid);
        nc_def_var(ncid, "lon", NC_DOUBLE, 1, &dimIds[1], &lon_varid);
        nc_def_var(ncid, "hs_interp", NC_DOUBLE, 2, dimIds, &interp_varid);
        nc_enddef(ncid);
        // write the data
        nc_put_var_double(ncid, lat_varid, &E->wavesLat[0]);
        nc_put_var_double(ncid, lon_varid, &E->wavesLon[0]);
        nc_put_var_double(ncid, interp_varid, nn_interp);
        // close the file
        nc_close(ncid);
        //exit(1);
        */


    } // end of loop over Tp time levels


    // estimate wave setup on roms
    // setup is only calculates at the coastal points
    // to calculate setup, we need Hs, Tp and slope at coastal pointers

    // read in the slope data

    get_coastal_slope(E);

    // malloc room for the setup field
    // jNOTE: fix up the size of the time dimension here!
	E->setup_on_roms = malloc3d_double(E->nTimeWavesSubset, E->nLonRho, E->nLatRho);
	// malloc room for the time interpolated data
	E->setup_on_roms_time_interp = malloc3d_double(E->nTimeRoms, E->nLonRho, E->nLatRho);
    // for each coastal point calculate the wave setup


	// assign closest slope value to costline derived from the roms rho_mask
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
				nearest_index[i][j] = get_nearest_slope_index(E, this_lat, this_lon, E->slopeLat, E->slopeLon);
			}
			else{
				//printf("fill it: i = %d, j = %d\n",i,j);
				nearest_index[i][j] = -999;
			}
		}
	}


	for(t=0;t<E->nTimeWavesSubset;t++){
		printf("#### t = %d\n",t);
		for(i=0;i<E->nLonRho;i++){
			for(j=0;j<E->nLatRho;j++){

				if(E->coastline_mask[i][j] == 1.0){ // if this point is a coastal point
					/*
					// get longitude and latitude of the point
					this_time = t;
					this_lat = E->lat_rho[i][j];
					this_lon = E->lon_rho[i][j];
					E->setup_on_roms[t][i][j] = get_nearest_setup(E, this_time, this_lat, this_lon, E->setup, E->wavesLat, E->wavesLon);
					//exit(1);
					*/

					E->setup_on_roms[t][i][j] = get_setup(E->Hs_on_roms[t][i][j], E->Tp_on_roms[t][i][j], E->slope[nearest_index[i][j]]);

				}
				else{
					//printf("fill it: i = %d, j = %d\n",i,j);
					E->setup_on_roms[t][i][j] = NC_FILL_DOUBLE;
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

    printf("creating %d interpolated time levels for the setup field\n", E->nTimeRoms);
	for(t=0;t<E->nTimeRoms;t++){
		// assume we are going to time interpolate
		interpolate = TRUE;
		time1 = -1;
		time2 = -1;

		// find the waveTimes that bound the current romsTime
		// get lower bound
		for(i=0;i<E->nTimeWaves;i++){
			if(E->wavesTime[i] <= E->romsTime[t])	// check this
				time1 = i;
		}
		// get upper bound
		for(i=E->nTimeWaves-1;i>=0;i--){
			if(E->wavesTime[i] >= E->romsTime[t])	// check this
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
		mu = ((E->romsTime[t] - E->wavesTime[time1]) / (E->wavesTime[time2] - E->wavesTime[time1]));
		printf("lower bound (wave time) = %f\n", E->wavesTime[time1] );
		printf("upper bound (wave time) = %f\n",  E->wavesTime[time2]);
		printf("target roms time = %f\n", E->romsTime[t]);
		printf("mu = %f\n",mu);
		for(i=0;i<E->nLonRho;i++){
			for(j=0;j<E->nLatRho;j++){
				// interp the wave setup onto the roms time
				// y1 is the first vale at the wave time
				y1 = E->setup_on_roms[time1][i][j];
				// y2 is the second value at the wave time
				y2 = E->setup_on_roms[time2][i][j];


				if(interpolate == TRUE)
					E->setup_on_roms_time_interp[t][i][j] = LinearInterpolate( y1, y2, mu);
				else{
                    //printf("no interp: t = %d, i = %d, j = %d\n", t,i,j);
                    //printf("\tsetup = %f\n", E->setup_on_roms[t][i][j]);
					E->setup_on_roms_time_interp[t][i][j] = E->setup_on_roms[t][i][j];

                }
			}
		}
	}

    free(E->Hs);
    free(E->Tp);
    free(E->Hs_on_roms);
    free(E->Tp_on_roms);
    free(E->wavesLon);
    free(E->wavesLat);
    free(E->setup_on_roms);
    free(E->nn_interp);
    free(E->nn_diff);
    free(E->slope);
    free(E->slopeLat);
    free(E->slopeLon);
}


double get_setup(double Hs, double Tp, double slope){

    double  L;
    double  setup;

    if(Hs == 0.0) Hs = 0.0;

    L = 9.8*pow(Tp,2.0)/(2.0*M_PI);
    setup = (0.31)*(0.8)*Hs*(0.32*pow(slope,1.0/7.0))*pow((Hs/L),-0.25);

    if(setup > 100.0){
        printf("WTF?? setup = %f, Hs = %f, Tp = %f, slope = %f\n", setup, Hs, Tp, slope);
    }
    return (setup);
}
