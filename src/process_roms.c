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
    double **cangle;

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

    //printf("ocean_time = %zu\n", E->nTimeRoms);

    // get the lat dimension sizes
    if((retval = nc_inq_dimid(ncid, "xi_rho", &varid)))
        fail("failed to get roms lat dimid: error is %d\n",retval);

    if((retval = nc_inq_dimlen(ncid,varid,&E->nLatRho)))
        fail("failed to get roms lat dimid: error is %d\n",retval);

    //printf("xi_rho = %zu\n", E->nLatRho);

    // get the lon dimension sizes
    if((retval = nc_inq_dimid(ncid, "eta_rho", &varid)))
        fail("failed to get roms lon_rho dimid: error is %d\n",retval);

    if((retval = nc_inq_dimlen(ncid,varid,&E->nLonRho)))
        fail("failed to read roms lon_rho data: error is %d\n",retval);

    //printf("eta_rho = %zu\n", E->nLonRho);

    // malloc room for the arrays
	E->romsTime = malloc(E->nTimeRoms*sizeof(double));
	E->lat_rho = malloc2d_double(E->nLonRho, E->nLatRho);
	E->lon_rho = malloc2d_double(E->nLonRho, E->nLatRho);
	E->mask_rho = malloc2d_double(E->nLonRho, E->nLatRho);
	E->zeta = malloc3d_double(E->nTimeRoms, E->nLonRho, E->nLatRho);
	E->bathymetry = malloc2d_double(E->nLonRho, E->nLatRho);
	E->db_dx = malloc2d_double(E->nLonRho, E->nLatRho);
	E->db_dy = malloc2d_double(E->nLonRho, E->nLatRho);
	E->bathymetry_gradient = malloc2d_double(E->nLonRho, E->nLatRho);
	E->pm = malloc2d_double(E->nLonRho, E->nLatRho);
	E->pn = malloc2d_double(E->nLonRho, E->nLatRho);
	E->bathymetry_slope = malloc2d_double(E->nLonRho, E->nLatRho);

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
	//printf("units = %s\n", E->roms_time_units);

	//printf("romsTime[0] = %f\n", E->romsTime[0]);
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
    //printf( "this date is %04d-%02d-%02d %02d:%02d:%06.3lf\n",yr, mo, day, hr, min, sec );



    nc_inq_varid(ncid, "lat_rho", &varid);
    if((retval = nc_get_var_double(ncid, varid, &E->lat_rho[0][0])))
		fail("failed to read roms lat_rho data: error is %d\n", retval);

		//printf("lat_rho[0][0] = %f\n", E->lat_rho[0][0]);

    nc_inq_varid(ncid, "lon_rho", &varid);
    if((retval = nc_get_var_double(ncid, varid, &E->lon_rho[0][0])))
		fail("failed to read roms lat_rho data: error is %d\n", retval);

	//printf("lon_rho[0][0] = %f\n", E->lon_rho[0][0]);

	// get the rho_mask
	nc_inq_varid(ncid, "mask_rho", &varid);
    if((retval = nc_get_var_double(ncid, varid, &E->mask_rho[0][0])))
		fail("failed to read roms mask_rho data: error is %d\n", retval);


	/*
	// get the bathymetry
	nc_inq_varid(ncid, "h", &varid);
    if((retval = nc_get_var_double(ncid, varid, &E->bathymetry[0][0])))
                fail("failed to read roms bathymetry data: error is %d\n", retval);


	// get pm
	nc_inq_varid(ncid, "pm", &varid);
    if((retval = nc_get_var_double(ncid, varid, &E->pm[0][0])))
                fail("failed to read roms pm data: error is %d\n", retval);

	// get pn
	nc_inq_varid(ncid, "pn", &varid);
    if((retval = nc_get_var_double(ncid, varid, &E->pn[0][0])))
                fail("failed to read roms pn data: error is %d\n", retval);
	*/

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

    //printf("spatial bounds from roms file:\n");
    //printf("\tlon: %f to %f\n", E->roms_min_lon, E->roms_max_lon);
    //printf("\tlat: %f to %f\n", E->roms_min_lat, E->roms_max_lat);



    // extract the coastline wet cells from the roms input
    E->coastline_mask = malloc2d_double(E->nLonRho, E->nLatRho);
	padded_mask = malloc2d_double(E->nLonRho+2, E->nLatRho+2);
	// initialize the field to be fill_value everywhere
	//printf("i size nLatRho = %d\n", E->nLatRho);
	//printf("j size nLonRho = %d\n", E->nLonRho);

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


   // open the ROMS grid file to read in the bathymetry
   // we need the bathymetry from this file so we can
   // calculate the gradient for wave setup

   if((retval = nc_open(E->grid_input, NC_NOWRITE, &ncid)))
        fail("failed to open roms grid file: error is %d\n",retval);

  // get the bathymetry
        nc_inq_varid(ncid, "bathymetry", &varid);
    if((retval = nc_get_var_double(ncid, varid, &E->bathymetry[0][0])))
                fail("failed to read roms bathymetry data: error is %d\n", retval);


        // get pm
        nc_inq_varid(ncid, "pm", &varid);
    if((retval = nc_get_var_double(ncid, varid, &E->pm[0][0])))
                fail("failed to read roms pm data: error is %d\n", retval);

        // get pn
        nc_inq_varid(ncid, "pn", &varid);
    if((retval = nc_get_var_double(ncid, varid, &E->pn[0][0])))
                fail("failed to read roms pn data: error is %d\n", retval); 


	nc_close(ncid);


	// process bathymetry data to make land cells zero
	// so we can calculate the gradient later
/*	
	// first find out what the 'depth' is for land masked points
	double min_depth = 5.0;
	for(i=0;i<E->nLonRho;i++){
                for(j=0;j<E->nLatRho;j++){
                        if(E->mask_rho[i][j] == 0.0){//assuming a min depth of 5 m
                                min_depth = E->bathymetry[i][j];
				break;
                        }
                }
        }
*/

	/*
	for(i=0;i<E->nLonRho;i++){
                for(j=0;j<E->nLatRho;j++){	
			if(E->bathymetry[i][j] < 0.0){//assuming a min depth of 5 m
				E->bathymetry[i][j] = 0.0;
			}
		}
	}
	//printf("bathy min depth = %f\n", min_depth);
	*/
	// now calculate bathymetry gradient
	// cal db/dx and db/dy using finite differences
	// the slope is then tan-1(magnutide(gradient))
	// slope angle = tan-1( sqrt(pow(db/dx,2)+pow(db/dy,2)) )
	//
	// backward difference for interior points

	// dh/dx 
	for(i=1;i<E->nLonRho;i++){
                for(j=0;j<E->nLatRho;j++){
			E->db_dx[i][j] = (E->bathymetry[i][j] - E->bathymetry[i-1][j])*E->pm[i][j];
		}
	}
	// get edge value using backward difference
	for(j=0;j<E->nLatRho;j++)
		E->db_dx[0][j] = (E->bathymetry[1][j] - E->bathymetry[0][j]) * E->pm[0][j];


	// dh/dy - backward difference
	for(i=0;i<E->nLonRho;i++){
                for(j=1;j<E->nLatRho;j++){
                        E->db_dy[i][j] = (E->bathymetry[i][j] - E->bathymetry[i][j-1])*E->pn[i][j];
                }
        }
	// get edge value using single sided difference
	for(i=0;i<E->nLonRho;i++)
		E->db_dy[i][0] = (E->bathymetry[i][1]-E->bathymetry[i][0])*E->pn[i][0];	


	for(i=0;i<E->nLonRho;i++){
                for(j=0;j<E->nLatRho;j++){
			E->bathymetry_gradient[i][j] = sqrt(pow(E->db_dx[i][j],2.0)+pow(E->db_dy[i][j],2.0));
			E->bathymetry_slope[i][j] = atan(E->bathymetry_gradient[i][j]);

                }
        }


  for(i=0;i<E->nLonRho;i++){
     for(j=0;j<E->nLatRho;j++){
        if(E->bathymetry_slope[i][j] < 0.001) 
              E->bathymetry_slope[i][j] = 0.001;
     }
  }


    // write gradient data to file
    // create the file
    int lon_varid, lat_varid;
    int lat_dimid, lon_dimid;
    int dimIds[2];
    int grad_varid, slope_varid, cangle_varid, dbdx_varid, dbdy_varid, h_varid;

    cangle = malloc2d_double(E->nLonRho, E->nLatRho);
    for(i=0;i<E->nLonRho;i++){
         for(j=0;j<E->nLatRho;j++){    
           if(E->coastline_mask[i][j] == 1)
		cangle[i][j] = E->bathymetry_slope[i][j];
           else
                cangle[i][j] = NC_FILL_DOUBLE;


	  if(E->mask_rho[i][j] == 0){
		E->bathymetry_gradient[i][j] = NC_FILL_DOUBLE;
		E->bathymetry_slope[i][j] = NC_FILL_DOUBLE; 
		E->bathymetry[i][j] = NC_FILL_DOUBLE;
		E->db_dx[i][j] = NC_FILL_DOUBLE;
		E->db_dy[i][j] = NC_FILL_DOUBLE;
	  }

         }
    }

    nc_create("bathy_gradient.nc", NC_CLOBBER, &ncid);
    // def dimensions
    nc_def_dim(ncid, "lat", E->nLatRho, &lat_dimid);
    nc_def_dim(ncid, "lon", E->nLonRho, &lon_dimid);
    // def vars
    dimIds[0] = lon_dimid;
    dimIds[1] = lat_dimid;
    nc_def_var(ncid, "lat", NC_DOUBLE, 2, dimIds, &lat_varid);
    nc_def_var(ncid, "lon", NC_DOUBLE, 2, dimIds, &lon_varid);
    nc_def_var(ncid, "gradient", NC_DOUBLE, 2, dimIds, &grad_varid);
    nc_def_var(ncid, "angle", NC_DOUBLE, 2, dimIds, &slope_varid);
    nc_def_var(ncid, "cangle", NC_DOUBLE, 2, dimIds, &cangle_varid);
    nc_def_var(ncid, "dbdx", NC_DOUBLE, 2, dimIds, &dbdx_varid);
    nc_def_var(ncid, "dbdy", NC_DOUBLE, 2, dimIds, &dbdy_varid);
    nc_def_var(ncid,"h", NC_DOUBLE,2,dimIds, &h_varid);
    nc_enddef(ncid);
    // write the data
    nc_put_var_double(ncid, lat_varid, &E->lat_rho[0][0]);
    nc_put_var_double(ncid, lon_varid, &E->lon_rho[0][0]);
    nc_put_var_double(ncid, h_varid, &E->bathymetry[0][0]);
    nc_put_var_double(ncid, grad_varid, &E->bathymetry_gradient[0][0]);
    nc_put_var_double(ncid, slope_varid, &E->bathymetry_slope[0][0]);
    nc_put_var_double(ncid, cangle_varid, &cangle[0][0]);
    nc_put_var_double(ncid, dbdx_varid, &E->db_dx[0][0]);
    nc_put_var_double(ncid, dbdy_varid, &E->db_dy[0][0]);
    // close the file
    nc_close(ncid);	

	


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

	free(E->zeta);

}
