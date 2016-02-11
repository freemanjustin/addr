#include "grid.h"

void process_roms(e *E){

    int	ncid;
	int varid;
	int retval;
	size_t attlen = 0;

    // for time conversion
	//static	char *calendar = "Standard";
    //ut_system	*u_system;
    double	sec;
    int	ierr, yr, mo, day, hr, min;


    int     i,j,t;
	int count;

    double **padded_mask;

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
    if( (E->roms_ref_time = ut_parse( E->u_system, E->roms_time_units, UT_ASCII )) == NULL ) {
        fprintf( stderr, "Error parsing units string \"%s\"\n", E->roms_time_units );
        exit(-1);
        }

    if( (ierr = utCalendar2_cal( E->romsTime[0], E->roms_ref_time, &yr, &mo, &day, &hr, &min, &sec, E->calendar )) != 0 ) {
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


    // find min and max longitude and latitudes for the roms data
    E->roms_min_lat = 9999.0; E->roms_max_lat = -99999.0;
    E->roms_min_lon = 9999.0; E->roms_max_lon = -99999.0;

    for(i=0;i<E->nLonRho;i++){
		for(j=0;j<E->nLatRho;j++){
			if(E->lon_rho[i][j] < E->roms_min_lon)
                E->roms_min_lon = E->lon_rho[i][j];
            if(E->lon_rho[i][j] > E->roms_max_lon)
                E->roms_max_lon = E->lon_rho[i][j];
            if(E->lat_rho[i][j] < E->roms_min_lat)
                E->roms_min_lat = E->lat_rho[i][j];
            if(E->lat_rho[i][j] > E->roms_max_lat)
                E->roms_max_lat = E->lat_rho[i][j];
		}
	}

    printf("spatial bounds from roms file:\n");
    printf("\tlon: %f to %f\n", E->roms_min_lon, E->roms_max_lon);
    printf("\tlat: %f to %f\n", E->roms_min_lat, E->roms_max_lat);



    // extract the coastline wet cells from the roms input
    E->coastline_mask = malloc2d_double(E->nLonRho, E->nLatRho);
	padded_mask = malloc2d_double(E->nLonRho+2, E->nLatRho+2);
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
	E->zeta_coast = malloc3d_double(E->nTimeRoms, E->nLonRho, E->nLatRho);
	for(t=0;t<E->nTimeRoms;t++){
		for(i=0;i<E->nLonRho;i++){
			for(j=0;j<E->nLatRho;j++){
                if(E->coastline_mask[i][j] == 1)
				    E->zeta_coast[t][i][j] = E->zeta[t][i][j];
                else
                    E->zeta_coast[t][i][j] = NC_FILL_DOUBLE;
			}
		}
	}

}
