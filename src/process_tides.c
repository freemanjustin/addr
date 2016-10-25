#include "grid.h"

void process_tides(e *E){

        int ncid;
        int varid;
        int retval;
        size_t attlen = 0;
        size_t from[4];
        size_t to[4];

        double sec;
        int ierr, yr, mo, day, hr, min;

        int i,j,t;
        int count;

        // read in the tide file
        if((retval = nc_open(E->tide_input, NC_NOWRITE, &ncid)))
                fail("failed to open tide input file: error is %d\n",retval);

        // get the lat dimension sizes
        if((retval = nc_inq_dimid(ncid, "lat", &varid)))
                fail("failed to get tide lat dimid: error is %d\n",retval);

        if((retval = nc_inq_dimlen(ncid,varid,&E->nLatTide)))
                fail("failed to get tide lat dim length: error is %d\n",retval);

        //printf("lat = %zu\n", E->nLatTide);

        // get the lon dimension sizes
        if((retval = nc_inq_dimid(ncid, "lon", &varid)))
                fail("failed to get tide lon dimid: error is %d\n",retval);

        if((retval = nc_inq_dimlen(ncid,varid,&E->nLonTide)))
                fail("failed to read tide lon dim length: error is %d\n",retval);

        //printf("tide lon = %zu\n", E->nLonTide);

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
        double dLatTide, dLonTide;
        dLatTide = fabs(E->tideLat[1] - E->tideLat[0])*2.0;
        dLonTide = fabs(E->tideLon[1] - E->tideLon[0])*2.0;
        for(i=0; i<E->nLatTide; i++) {
                if(E->tideLat[i] < (E->roms_min_lat-dLatTide))
                        lat_start = i;
        }
        int lat_end;
        for(i=E->nLatTide-1; i>=0; i--) {
                if(E->tideLat[i] > (E->roms_max_lat+dLatTide))
                        lat_end = i;
        }
        //printf("tide data start lat = %f (%d), end lat = %f (%d)\n", E->tideLat[lat_start],lat_start, E->tideLat[lat_end],lat_end);
        int lon_start;
        for(i=0; i<E->nLonTide; i++) {
                if(E->tideLon[i] < (E->roms_min_lon-dLonTide))
                        lon_start = i;
        }
        int lon_end;
        for(i=E->nLonTide-1; i>=0; i--) {
                if(E->tideLon[i] > (E->roms_max_lon+dLonTide))
                        lon_end = i;
        }

        //printf("tide data start lon = %f, end lon = %f\n", E->tideLon[lon_start], E->tideLon[lon_end]);

        // TODO: add some error checking to the bounds code.
        // for example, if the spatial extent does not overlap then throw an exception

        // now just read in what we want from the files
        free(E->tideLat);
        free(E->tideLon);

	// gap testing
	lat_start = 0;
	lat_end = E->nLatTide-1;
	lon_start = 0;
	lon_end = E->nLonTide-1;

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

        //printf("tide lev = %zu\n", E->nLevTide);

        // get the time dimension sizes
        if((retval = nc_inq_dimid(ncid, "time", &varid)))
                fail("failed to get tide tide dimid: error is %d\n",retval);

        if((retval = nc_inq_dimlen(ncid,varid,&E->nTimeTide)))
                fail("failed to read tide time data: error is %d\n",retval);

        //printf("tide time = %zu\n", E->nTimeTide);

        // malloc room for the arrays
        //(time, lev, lat, lon)
        E->tideTime = malloc(E->nTimeTide*sizeof(double));

        // read in the time data
        nc_inq_varid(ncid, "time", &varid);
        if((retval = nc_get_var_double(ncid, varid, &E->tideTime[0])))
                fail("failed to read tide time data data: error is %d\n",retval);
        // get the time metadata units
        nc_inq_attlen (ncid, varid, "units", &attlen);
        E->tide_time_units = (char *) malloc(attlen + 1); /* + 1 for trailing null */
        nc_get_att_text(ncid, varid, "units", E->tide_time_units);
        E->tide_time_units[attlen] = '\x0';
        //printf("tide time units = %s\n", E->tide_time_units);

        // Make the Calendar calls
        // Parse the units strings
        if( (E->tide_ref_time = ut_parse( E->u_system, E->tide_time_units, UT_ASCII )) == NULL ) {
                fprintf( stderr, "Error parsing units string \"%s\"\n", E->roms_time_units );
                exit(-1);
        }

        // put the tide time on the same time units as the roms time
        for(t=0; t<E->nTimeTide; t++) {

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
        for(t=0; t<E->nTimeTide; t++) {
                if(E->tideTime[t]<=E->start_time_roms)
                        E->tide_start_time_index = t;
        }
        if(E->tide_start_time_index == -1) {
                fprintf(stderr,"couldn't find a matching start time in the tide file.\n");
                fprintf(stderr,"check to make sure the tide file times sufficiently overlap\n");
                fprintf(stderr,"the ROMS times.\n\n");
                exit(1);
        }

        // get the end index for the tide file
        E->tide_end_time_index = -1;
        E->end_time_roms = E->romsTime[E->nTimeRoms-1];
        for(t=E->nTimeTide-1; t>=0; t--) {
                //printf("t = %d, tide_time = %f\n",t,E->tideTime[t]);
                if(E->tideTime[t] >= E->end_time_roms)
                        E->tide_end_time_index = t;
        }

        if(E->tide_end_time_index == -1) {
                fprintf(stderr,"couldn't find a matching end time in the tide file.\n");
                fprintf(stderr,"check to make sure the tide file times sufficiently overlap\n");
                fprintf(stderr,"the ROMS times.\n\n");
                fprintf(stderr,"end time ROMS = %f\n", E->end_time_roms);
                fprintf(stderr,"end time tide = %f\n", E->tideTime[E->nTimeTide-1]);
                exit(1);
        }

        //printf("Tide start index = %d\n", E->tide_start_time_index);
        //printf("Tide end index = %d\n", E->tide_end_time_index);

        E->nTimeTideSubset = (E->tide_end_time_index-E->tide_start_time_index+1);

        // make the time vector for the output file
        E->tide_interp_time = malloc(E->nTimeTideSubset*sizeof(double));
        count=0;
        for(t=E->tide_start_time_index; t<E->tide_end_time_index; t++) {
                E->tide_interp_time[count] = E->tideTime[t];
                count++;
        }


        // read the data
        E->tide_data = malloc4d_double(E->nTimeTideSubset, E->nLevTide, E->nLatTide, E->nLonTide);

        from[0] = E->tide_start_time_index;     to[0] = E->nTimeTideSubset;
        from[1] = 0;                            to[1] = E->nLevTide;
        from[2] = lat_start;                    to[2] = lat_end - lat_start;//E->nLatTide;
        from[3] = lon_start;                    to[3] = lon_end - lon_start;//E->nLonTide;

        nc_inq_varid(ncid, "z", &varid);
        // should only read in the data we need
        if((retval = nc_get_vara_double(ncid, varid, from, to, &E->tide_data[0][0][0][0])))
                fail("failed to read tide data data: error is %d\n",retval);

        // close the tide file
        nc_close(ncid);

    #ifdef CHECK
        // temporarily mess with this input data to check!
        for(t=0; t<E->nTimeTideSubset; t++) {
                for(i=0; i<E->nLatTide; i++) {
                        for(j=0; j<E->nLonTide; j++) {
                                E->tide_data[t][0][i][j] = (double)t;
                        }
                }
        }
    #endif

        // do a natural neighbour interpolation on the tide data we just read in
        E->tide_on_roms = malloc3d_double(E->nTimeTideSubset, E->nLonRho, E->nLatRho);

        // nn optimized working arrays
        E->pin = malloc(E->nLonTide * E->nLatTide * sizeof(point));
        E->zin = malloc(E->nLonTide * E->nLatTide * sizeof(double));

        E->xout = malloc(E->nLonRho * E->nLatRho * sizeof(double));
        E->yout = malloc(E->nLonRho * E->nLatRho * sizeof(double));
        E->zout = malloc(E->nLonRho * E->nLatRho * sizeof(double));


        // find out how many valid data points we have
        // and setup the input array for lib-nn
        //time_vector[t] = (double)t;
        //printf("setting up source grid for nn...\n");
        E->nin = 0;
        for(i=0; i<E->nLatTide; i++) {
                for(j=0; j<E->nLonTide; j++) {
                        if(E->tide_data[0][0][i][j] > -999.0) {
                                E->pin[E->nin].x = E->tideLon[j];
                                E->pin[E->nin].y = E->tideLat[i];
                                E->nin++;
                        }
                }
        }
        //printf("done\n");fflush(stdout);


        // now set up the output array for the nn interpolation
        // this is the roms grid
        E->nout = 0;
        for(i=0; i<E->nLonRho; i++) {
                for(j=0; j<E->nLatRho; j++) {
                      E->xout[E->nout] = E->lon_rho[i][j];
                      E->yout[E->nout] = E->lat_rho[i][j];
                      E->zout[E->nout] = NaN;
                      E->nout++;
                }
        }
        //printf("done\n");fflush(stdout);

        // setup the natural neighbour interpolation
        // only need to do this once
        E->d = delaunay_build(E->nin, E->pin, 0, NULL, 0, NULL);

        // create interpolator
        E->nn = nnai_build(E->d, E->nout, E->xout, E->yout);

        // for each time level
        for(t=0; t<E->nTimeTideSubset; t++) {

                // setup nodal values for the nn interps
                E->nin = 0;
                for(i=0; i<E->nLatTide; i++) {
                        for(j=0; j<E->nLonTide; j++) {
                                if(E->tide_data[t][0][i][j] > -999.0) {
                                        E->pin[E->nin].z = E->tide_data[t][0][i][j];
                                        E->zin[E->nin] = E->pin[E->nin].z;
                                        E->nin++;
                                }
                        }
                }

                // do the interpolation
                nnai_interpolate(E->nn, E->zin, E->zout);

                // splat interpolated values onto the roms grid and apply land-sea mask
                count = 0;
                for(i=0; i<E->nLonRho; i++) {
                        for(j=0; j<E->nLatRho; j++) {
                                if(E->coastline_mask[i][j] == 0){
                                //if(E->mask_rho[i][j] == 0){
                                        E->tide_on_roms[t][i][j] = NC_FILL_DOUBLE;
                                }
                                else{
                                        E->tide_on_roms[t][i][j] = E->zout[count];
                                }
                                count++;
                        }
                }

                /*
                // write it out to check
                //E->nn_nx * E->nn_ny, E->nn_interp,
                int lat_dimid, lon_dimid, time_dimid, dimIds[2];
                int lat_varid, lon_varid, time_varid, interp_varid;
                // create the file
                nc_create("tide_interp.nc", NC_CLOBBER, &ncid);
                // def dimensions
                nc_def_dim(ncid, "lat", E->nLonRho, &lat_dimid);
                nc_def_dim(ncid, "lon", E->nLatRho, &lon_dimid);
                // def vars
                dimIds[0] = lat_dimid;
                dimIds[1] = lon_dimid;
                //nc_def_var(ncid, "lat", NC_DOUBLE, 1, &dimIds[0], &lat_varid);
                //nc_def_var(ncid, "lon", NC_DOUBLE, 1, &dimIds[1], &lon_varid);
                nc_def_var(ncid, "tide_interp_on_roms", NC_DOUBLE, 2, dimIds, &interp_varid);
                nc_enddef(ncid);
                // write the data
                //nc_put_var_double(ncid, lat_varid, &E->wavesLat[0]);
                //nc_put_var_double(ncid, lon_varid, &E->wavesLon[0]);
                nc_put_var_double(ncid, interp_varid, &E->tide_on_roms[0][0][0]);
                // close the file
                nc_close(ncid);
                exit(1);
                */

        } // end of loop over tide time levels
        free(E->pin);
        free(E->zin);

        free(E->xout);
        free(E->yout);
        free(E->zout);

        free(E->d);
        free(E->nn);

        // time interp tide data to roms time

        // time interpolate the wavesetup data onto the roms time vector
        //printf("creating %d interpolated time levels for the tide field\n", E->nTimeRoms);
        E->tide_on_roms_time_interp = malloc3d_double(E->nTimeRoms, E->nLonRho, E->nLatRho);
        // initialize this array
        for(t=0; t<E->nTimeRoms; t++) {
                for(i=0; i<E->nLonRho; i++) {
                        for(j=0; j<E->nLatRho; j++) {
                                E->tide_on_roms_time_interp[t][i][j] = NC_FILL_DOUBLE;
                        }
                }
        }
        // target time vector = E->romsTime
        // source time vector = E->wavesTime
        // y value vector
        double *ypts = malloc(E->nTimeTideSubset*sizeof(double));
        double *interp_y = malloc(E->nTimeRoms*sizeof(double));
        for(i=0; i<E->nLonRho; i++) {
                for(j=0; j<E->nLatRho; j++) {
                        if(E->tide_on_roms[0][i][j] != NC_FILL_DOUBLE) {
                                for(t=0; t<E->nTimeTideSubset; t++) {
                                        // get the wave setup vector at this location
                                        ypts[t] = E->tide_on_roms[t][i][j];
                                }
                                time_interp_field(&E->tideTime[E->tide_start_time_index], &ypts[0], E->nTimeTideSubset, &E->romsTime[0], &interp_y[0], E->nTimeRoms);
                                // now put this data into the time interp array
                                for(t=0; t<E->nTimeRoms; t++) {
                                        // get the wave setup vector at this location
                                        E->tide_on_roms_time_interp[t][i][j] = interp_y[t];
                                }
                        }
                }
        }
        free(ypts);
        free(interp_y);

        free(E->tideLon);
        free(E->tideLat);
        free(E->tide_data);

        free(E->tide_on_roms);
}
