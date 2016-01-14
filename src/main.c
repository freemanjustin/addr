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

	// for time conversion
	static	char *calendar = "Standard";
    ut_system	*u_system;
    ut_unit	*u_1;
    double	tval, tval_inv, sec;
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
	/* Initialize the udunits-2 library */
    ut_set_error_message_handler(ut_ignore);
    if( (u_system = ut_read_xml( NULL )) == NULL ) {
        fprintf( stderr, "Error initializing udunits-2 unit system\n" );
        exit(-1);
        }
    ut_set_error_message_handler(ut_write_to_stderr);

    // read input XML
    get_params(E);

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
	    if( (u_1 = ut_parse( u_system, E->roms_time_units, UT_ASCII )) == NULL ) {
	        fprintf( stderr, "Error parsing units string \"%s\"\n", E->roms_time_units );
	        exit(-1);
	        }

	    if( (ierr = utCalendar2_cal( E->romsTime[0], u_1, &yr, &mo, &day, &hr, &min, &sec, calendar )) != 0 ) {
	        fprintf( stderr, "Error on utCalendar2_cal call: %s\n", ccs_err_str(ierr) );
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
	E->tideLon = malloc(E->nLonTide*sizeof(double));
	E->tideLat = malloc(E->nLatTide*sizeof(double));
    E->tide_data = malloc4d_double(E->nTimeTide, E->nLevTide, E->nLatTide, E->nLonTide);


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
	// read in the number of lons and lon data values from the tide file

	// read in the number of lats and lat data values from the tide file

	// read in the tide data from the input file

	// interpolate the tide data for each lon_rho and lat_rho point

	// target grid for tide data
	// this needs a time dimension from the roms file
	E->tide_on_roms = malloc2d_double(E->nLatRho, E->nLonRho);
	for(i=0;i<E->nLatRho;i++)
		for(j=0;j<E->nLonRho;j++)
			E->tide_on_roms[i][j] = 0.0;
	// close the tide file
	nc_close(ncid);

	// interpolate tide to roms grid
    interp_tide_to_roms(E);
	// write the interpolated field to file
	write_netcdf(E);

	exit(1);

    // grid settings
    E->g.nX = E->g.X/E->g.resol + 1;
    E->g.nY = E->g.Y/E->g.resol + 1;

    E->g.rotangle = E->g.rotangle/180.0*pi; // Convert Angle for grid rotation from degrees to radians
    latdist  = spheriq_dist(E->g.lon,E->g.lat,E->g.lon,E->g.lat+1.0, 0); // Length (in meters) of 1 degree of latitude

    // malloc arrays
    malloc_arrays(E);

    E->x[0] = 0.0;
    for(i=1;i<E->g.nX;i++){
        E->x[i] = E->x[i-1]+E->g.resol;
    }
    E->y[0] = 0.0;
    for(i=1;i<=E->g.nY;i++){
        E->y[i] = E->y[i-1]+E->g.resol;
    }

    Lm = E->g.nX-2;
    Mm = E->g.nY-2;
    Lp = Lm + 2;
    Mp = Mm + 2;
    L = Lp - 1;
    M = Mp -1;


    // RHO GRID
    // Create non-georeferenced grid in meters (origin = 0,0)
    for(i=0;i<=E->g.nY;i++){
        for(j=0;j<E->g.nX;j++){
            E->x_rho[i][j] = E->x[j];
            E->y_rho[i][j] = E->y[i];
        }
    }

    // Rotate grid
    for(i=0;i<=E->g.nY;i++){
        for(j=0;j<E->g.nX;j++){
            E->Rx_rho[i][j] = E->x_rho[i][j] * cos(E->g.rotangle) - E->y_rho[i][j] * sin(E->g.rotangle);
            E->Ry_rho[i][j] = E->x_rho[i][j] * sin(E->g.rotangle) + E->y_rho[i][j] * cos(E->g.rotangle);
        }
    }

    // Estimate Latitude and Longitude of each Grid point
    E->lat_rho[0][0] = E->g.lat + ( E->Ry_rho[0][0] / latdist);
    E->lon_rho[0][0] = E->g.lon + (E->Rx_rho[0][0]/spheriq_dist(E->g.lon,E->lat_rho[0][0],E->g.lon+1,E->lat_rho[0][0], 0));

        for(i=1;i<=E->g.nY;i++){
        E->lat_rho[i][0] = E->g.lat + (E->Ry_rho[i][0]/ latdist);
        E->lon_rho[i][0] = E->lon_rho[i-1][0] + ((E->Rx_rho[i][0]-E->Rx_rho[i-1][0]) / spheriq_dist(E->lon_rho[i-1][0],E->lat_rho[i][0],E->lon_rho[i-1][0]+1,E->lat_rho[i][0], 0));
    }

    for(i=0;i<=E->g.nY;i++){
        for(j=1;j<E->g.nX;j++){
            E->lat_rho[i][j] = E->g.lat + (E->Ry_rho[i][j] / latdist);
            E->lon_rho[i][j] = E->lon_rho[i][j-1] + ((E->Rx_rho[i][j]-E->Rx_rho[i][j-1]) / spheriq_dist(E->lon_rho[i][j-1],E->lat_rho[i][j],E->lon_rho[i][j-1]+1,E->lat_rho[i][j],0));
        }
    }

    free(E->Rx_rho);
    free(E->Ry_rho);



    // U GRID
    // Create non-georeferenced grid in meters (origin = 0,0)

    for(i=0;i<E->g.nY;i++){
        for(j=0;j<E->g.nX-1;j++){
            E->x_u[i][j] = 0.5*(E->x_rho[i][j]+E->x_rho[i][j+1]);
            E->y_u[i][j] = 0.5*(E->y_rho[i][j]+E->y_rho[i][j+1]);
        }
    }


    // Rotate grid
    for(i=0;i<E->g.nY;i++){
        for(j=0;j<E->g.nX-1;j++){
            E->Rx_u[i][j] = E->x_u[i][j] * cos(E->g.rotangle) - E->y_u[i][j] * sin(E->g.rotangle);
            E->Ry_u[i][j] = E->x_u[i][j] * sin(E->g.rotangle) + E->y_u[i][j] * cos(E->g.rotangle);

        }
    }

    // Estimate Latitude and Longitude of each Grid point
    E->lat_u[0][0] = E->g.lat + (E->Ry_u[0][0] / latdist);
    E->lon_u[0][0] = E->g.lon + (E->Rx_u[0][0] / spheriq_dist(E->g.lon,E->lat_u[0][0],E->g.lon+1,E->lat_u[0][0],0));

    for(i=1;i<E->g.nY;i++){
        E->lat_u[i][0] = E->g.lat + (E->Ry_u[i][0]/ latdist);
        E->lon_u[i][0] = E->lon_u[i-1][0] + ((E->Rx_u[i][0]-E->Rx_u[i-1][0]) / spheriq_dist(E->lon_u[i-1][0],E->lat_u[i][0],E->lon_u[i-1][0]+1,E->lat_u[i][0],0));
    }

    for(i=0;i<E->g.nY;i++){
        for(j=1;j<E->g.nX-1;j++){
            E->lat_u[i][j] = E->g.lat + (E->Ry_u[i][j] / latdist);
            E->lon_u[i][j] = E->lon_u[i][j-1] + ((E->Rx_u[i][j]-E->Rx_u[i][j-1]) / spheriq_dist(E->lon_u[i][j-1],E->lat_u[i][j],E->lon_u[i][j-1]+1,E->lat_u[i][j],0));
        }
    }

    free(E->Rx_u);
    free(E->Ry_u);
    free(E->x_u);
    free(E->y_u);


    // Calculate angle variable - angle at rho points
    for(i=0;i<E->g.nY;i++){
          for(j=0;j<E->g.nX;j++){

            E->angle[i][j] = bearing(E->lat_rho[i][j],E->lon_rho[i+1][j],E->lat_rho[i+1][j], E->lon_rho[i][j] );

        }
    }

    // V GRID
    // Create non-georeferenced grid in meters (origin = 0,0)

    for(i=0;i<E->g.nY-1;i++){
        for(j=0;j<E->g.nX;j++){
            E->x_v[i][j] = 0.5*(E->x_rho[i][j]+E->x_rho[i+1][j]);
            E->y_v[i][j] = 0.5*(E->y_rho[i][j]+E->y_rho[i+1][j]);
        }
    }

    // Rotate grid

    for(i=0;i<E->g.nY-1;i++){
        for(j=0;j<E->g.nX;j++){
            E->Rx_v[i][j] = E->x_v[i][j] * cos(E->g.rotangle) - E->y_v[i][j] * sin(E->g.rotangle);
            E->Ry_v[i][j] = E->x_v[i][j] * sin(E->g.rotangle) + E->y_v[i][j] * cos(E->g.rotangle);

        }
    }

    // Estimate Latitude and Longitude of each Grid point
    E->lat_v[0][0] = E->g.lat + (E->Ry_v[0][0] / latdist);
    E->lon_v[0][0] = E->g.lon + (E->Rx_v[0][0] / spheriq_dist(E->g.lon,E->lat_v[0][0],E->g.lon+1,E->lat_v[0][0],0));

    for(i=1;i<E->g.nY-1;i++){
        E->lat_v[i][0] = E->g.lat + (E->Ry_v[i][0]/ latdist);
        E->lon_v[i][0] = E->lon_v[i-1][0] + ((E->Rx_v[i][0]-E->Rx_v[i-1][0]) / spheriq_dist(E->lon_v[i-1][0],E->lat_v[i][0],E->lon_v[i-1][0]+1,E->lat_v[i][0],0));
    }


    for(i=0;i<E->g.nY-1;i++){
        for(j=1;j<E->g.nX;j++){
            E->lat_v[i][j] = E->g.lat + (E->Ry_v[i][j] / latdist);
            E->lon_v[i][j] = E->lon_v[i][j-1] + ((E->Rx_v[i][j]-E->Rx_v[i][j-1]) / spheriq_dist(E->lon_v[i][j-1],E->lat_v[i][j],E->lon_v[i][j-1]+1,E->lat_v[i][j],0));
        }
    }

    free(E->Rx_v);
    free(E->Ry_v);
    free(E->x_v);
    free(E->y_v);

    // PSI GRID
    // Create non-georeferenced grid in meters (origin = 0,0)

    for(i=0;i<E->g.nY-1;i++){
        for(j=0;j<E->g.nX-1;j++){
            E->x_psi[i][j] = 0.5*(E->x_rho[i][j]+E->x_rho[i+1][j+1]);
            E->y_psi[i][j] = 0.5*(E->y_rho[i][j]+E->y_rho[i+1][j+1]);
        }
    }

    // Rotate grid

    for(i=0;i<E->g.nY-1;i++){
        for(j=0;j<E->g.nX-1;j++){
            E->Rx_psi[i][j] = E->x_psi[i][j] * cos(E->g.rotangle) - E->y_psi[i][j] * sin(E->g.rotangle);
            E->Ry_psi[i][j] = E->x_psi[i][j] * sin(E->g.rotangle) + E->y_psi[i][j] * cos(E->g.rotangle);

        }
    }

    // Estimate Latitude and Longitude of each Grid point

    E->lat_psi[0][0] = E->g.lat + (E->Ry_psi[0][0] / latdist);
    E->lon_psi[0][0] = E->g.lon + (E->Rx_psi[0][0] / spheriq_dist(E->g.lon,E->lat_psi[0][0],E->g.lon+1,E->lat_psi[0][0],0));

    for(i=1;i<E->g.nY-1;i++){
        E->lat_psi[i][0] = E->g.lat + (E->Ry_psi[i][0]/ latdist);
        E->lon_psi[i][0] = E->lon_psi[i-1][0] + ((E->Rx_psi[i][0]-E->Rx_psi[i-1][0]) / spheriq_dist(E->lon_psi[i-1][0],E->lat_psi[i][0],E->lon_psi[i-1][0]+1,E->lat_psi[i][0],0));
    }


    for(i=0;i<E->g.nY-1;i++){
        for(j=1;j<E->g.nX-1;j++){
            E->lat_psi[i][j] = E->g.lat + (E->Ry_psi[i][j] / latdist);
            E->lon_psi[i][j] = E->lon_psi[i][j-1] + ((E->Rx_psi[i][j]-E->Rx_psi[i][j-1]) / spheriq_dist(E->lon_psi[i][j-1],E->lat_psi[i][j],E->lon_psi[i][j-1]+1,E->lat_psi[i][j],0));
        }
    }

    free(E->x_rho);
    free(E->y_rho);
    free(E->x_psi);
    free(E->y_psi);
    free(E->Rx_psi);
    free(E->Ry_psi);


    // Grid spacing and other grid parameters

    el = E->lat_u[E->g.nY-1][0] - E->lat_u[0][0];
    xl = E->lon_v[0][E->g.nX-1] - E->lon_v[0][0];

    for(i=0;i<E->g.nY;i++){
        for(j=1;j<E->g.nX-1;j++){
            E->dx[i][j] = spheriq_dist(E->lon_u[i][j],E->lat_u[i][j],E->lon_u[i][j-1],E->lat_u[i][j-1],0);
        }
    }

    for(i=0;i<E->g.nY;i++){
        E->dx[i][0] = E->dx[i][1];
        E->dx[i][E->g.nX-1] = E->dx[i][E->g.nX-2];
    }

    for(i=1;i<E->g.nY-1;i++){
        for(j=0;j<E->g.nX;j++){
            E->dy[i][j] = spheriq_dist(E->lon_v[i][j],E->lat_v[i][j],E->lon_v[i-1][j],E->lat_v[i-1][j],0);
        }
    }

    for(i=0;i<E->g.nX;i++){
        E->dy[0][i] = E->dy[1][i];
        E->dy[E->g.nY-1][i] = E->dy[E->g.nY-2][i];
    }

    for(i=0;i<E->g.nY;i++){
        for(j=0;j<E->g.nX;j++){
            E->pm[i][j] = 1.0/E->dx[i][j];
            //printf("pm[%d][%d] = %g\n", i, j, pm[i][j]);
            E->pn[i][j] = 1.0/E->dy[i][j];
        }
    }

    for(i=0;i<E->g.nY;i++){
        for(j=0;j<E->g.nX;j++){
            E->dndx[i][j] = E->dmde[i][j] = 0.0;
        }
    }

    for(i=1;i<E->g.nY-1;i++){
        for(j=1;j<E->g.nX-1;j++){
            E->dndx[i][j] = 0.5*(1.0/E->pn[i][j+1] - 1.0/E->pn[i][j-1]);
            E->dmde[i][j] = 0.5*(1.0/E->pm[i+1][j] - 1.0/E->pm[i-1][j]);
        }
    }

    // Coriolis
    // f = 2 .* 7.29E-5 .* sin(lat_rho .* (pi/180)); %Estimation of Coriolis over the grid domain. OMEGA=7.29E-5
    // More info: http://en.wikipedia.org/wiki/Coriolis_effect#Formula
    for(i=0;i<E->g.nY;i++){
        for(j=0;j<E->g.nX;j++){
           E->f[i][j] = 2.0 * 7.29E-5 * sin(E->lat_rho[i][j] * (M_PI/180.0));
        }
    }

    // interpolate bathymetry
    interp_bathy_on_grid(E);
    // free bathymetry memory
    free(E->b.field);
    free(E->b.lon);
    free(E->b.lat);

    // generate land sea mask
    // 0 == land
    // 1 == water

    // Land/Sea mask on RHO-points
    // h[i][j] is now positive
    for(i=0;i<E->g.nY;i++){
        for(j=0;j<E->g.nX;j++){
            // original
            //E->mask_rho[i][j] = E->h[i][j] < E->b.min_depth ? 0 : 1;
            // psandery mod - make wet cells in the bathymetry wet in the mask
            E->mask_rho[i][j] = E->h[i][j] < 0.0 ? 0 : 1;
        }
    }

    // Land/Sea mask on U-points.
    for(i=0;i<E->g.nY;i++){
        for(j=0;j<E->g.nX-1;j++){
            E->mask_u[i][j] = E->mask_rho[i][j] * E->mask_rho[i][j+1];
        }
    }

    //  Land/Sea mask on V-points.
    for(i=0;i<E->g.nY-1;i++){
        for(j=0;j<E->g.nX;j++){
            E->mask_v[i][j] = E->mask_rho[i][j] * E->mask_rho[i+1][j];
        }
    }

    // Land/Sea mask on PSI-points.
    for(i=0;i<E->g.nY-1;i++){
        for(j=0;j<E->g.nX-1;j++){
            E->mask_psi[i][j] = E->mask_rho[i][j] * E->mask_rho[i+1][j];
        }
    }

    // apply min depth to bathymetry
    for(i=0;i<E->g.nY;i++){
        for(j=0;j<E->g.nX;j++){
            E->h[i][j] = E->h[i][j] < E->b.min_depth ? E->b.min_depth : E->h[i][j];
            // make bathymetry positive for ROMS
            //E->h[i][j] = fabs(E->h[i][j]);
        }
    }

    // smooth bathymetry?

    // save netcdf grid file
    write_netcdf(E);

	return 0;
}
