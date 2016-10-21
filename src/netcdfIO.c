// gridr
//
// freeman.justin@gmail.com



#include "grid.h"



void write_netcdf(e *E){

	//float	fillValue = -1e34;
	float fillValue = NC_FILL_FLOAT;


	create_netcdf(E, E->fname, &E->ncid);

    defdims_netcdf(E);

	// setup variables
	defvars(E);


	/*
	// add global metadata
	add_global_metadata(E, E->ncid);
	*/


    // write the data to the netcdf file
    write_data(E);

    nc_close(E->ncid);
}


void create_netcdf(e *E, char *fname, int *ncid){

	//int	old_fill_mode;

	if ( (E->retval = nc_create(fname, NC_CLOBBER, ncid) ) )
		fail("couldn't create netcdf out file %s\n",fname);

	// set the fill mode to be NO_FILL
	//nc_set_fill(*ncid, NC_NOFILL, &old_fill_mode);

}

void defdims_netcdf(e *E){


	E->one = 1;

	// define the dimensions
	if ((E->retval = nc_def_dim(E->ncid, "time", E->nTimeRoms, &E->ocean_time_dimid)))
		fail("nc_def_dim failed!\n");

	if ((E->retval = nc_def_dim(E->ncid, "xi_rho", E->nLatRho, &E->xi_rho_dimid)))
		fail("nc_def_dim failed!\n");

    if ((E->retval = nc_def_dim(E->ncid, "eta_rho", E->nLonRho, &E->eta_rho_dimid)))
		fail("nc_def_dim failed!\n");


	// end define mode for this file
	if ((E->retval = nc_enddef(E->ncid)))
      	fail("nc_enddef failed\n");

}


void defvars(e *E){

    // setup dimids


		E->dimIdsTide[0] = E->ocean_time_dimid;
		E->dimIdsTide[1] = E->eta_rho_dimid;
    E->dimIdsTide[2] = E->xi_rho_dimid;

    E->dimIdsRho[0] = E->eta_rho_dimid;
    E->dimIdsRho[1] = E->xi_rho_dimid;


	defvar_netcdf(E, E->ncid, "tide", NC_DOUBLE, 3, &E->dimIdsTide[0], &E->vid_tide);
    add_txt_attribute_netcdf(E, E->ncid, E->vid_h, "long_name", "tide at RHO-points");
		add_txt_attribute_netcdf(E, E->ncid, E->vid_h, "units", "m");




	defvar_netcdf(E, E->ncid, "lat", NC_DOUBLE, 2, E->dimIdsRho, &E->vid_lat_rho);
    add_txt_attribute_netcdf(E, E->ncid, E->vid_lat_rho, "long_name", "latitude");
		add_txt_attribute_netcdf(E, E->ncid, E->vid_lat_rho, "standard_name", "latitude");
		add_txt_attribute_netcdf(E, E->ncid, E->vid_lat_rho, "units", "degrees north");
    //add_txt_attribute_netcdf(E, E->ncid, E->vid_lat_rho, "_CoordinateAxisType", "Lat");
		add_txt_attribute_netcdf(E, E->ncid, E->vid_lat_rho, "axis", "Y");

	defvar_netcdf(E, E->ncid, "lon", NC_DOUBLE, 2, &E->dimIdsRho[0], &E->vid_lon_rho);
    add_txt_attribute_netcdf(E, E->ncid, E->vid_lon_rho, "long_name", "longitude");
		add_txt_attribute_netcdf(E, E->ncid, E->vid_lon_rho, "standard_name", "longitude");
		add_txt_attribute_netcdf(E, E->ncid, E->vid_lon_rho, "units", "degrees east");
    //add_txt_attribute_netcdf(E, E->ncid, E->vid_lon_rho, "_CoordinateAxisType", "Lon");
		add_txt_attribute_netcdf(E, E->ncid, E->vid_lon_rho, "axis", "X");


	defvar_netcdf(E, E->ncid, "time", NC_DOUBLE, 1, &E->dimIdsTide[0], &E->vid_ocean_time);
    add_txt_attribute_netcdf(E, E->ncid, E->vid_ocean_time, "long_name", "time");
		add_txt_attribute_netcdf(E, E->ncid, E->vid_ocean_time, "units", "seconds since 1970-01-01 00:00:00");
    add_txt_attribute_netcdf(E, E->ncid, E->vid_ocean_time, "calendar", "gregorian");
		add_txt_attribute_netcdf(E, E->ncid, E->vid_ocean_time, "_CoordinateAxisType", "time");
		add_txt_attribute_netcdf(E, E->ncid, E->vid_ocean_time, "axis", "T");



}



void defvar_netcdf(e *E, int ncid, char *var_name, nc_type type, int dims, int *dimIds, int *varid){


	if((E->retval = nc_redef(ncid)))
		fail("nc_redef failed\n");

	if ((E->retval = nc_def_var(ncid, var_name, type, dims, dimIds, varid)))
		fail("nc_def_var failed\n");

	if ((E->retval = nc_enddef(ncid)))
		fail("nc_enddef failed\n");

}


void add_txt_attribute_netcdf(e *E, int ncid, int varid, char* att_name, char* att_value){

	if((E->retval = nc_redef(ncid)))
		fail("nc_redef failed\n");

	if(( E->retval = nc_put_att_text(ncid, varid, att_name, strlen(att_value), att_value)))
		fail("nc_put_att_txt failed\n");

	if ((E->retval = nc_enddef(ncid)))
		fail("nc_enddef failed\n");
}


void add_global_metadata(e *E, int ncid){

    // TODO!
}


void write_data(e *E){


    // write tide
    if ((E->retval = nc_put_var_double(E->ncid, E->vid_tide, &E->tide_on_roms[0][0][0])))
        fail("put_var_ failed. Error code = %d\n",E->retval);


    // lat_rho
    if ((E->retval = nc_put_var_double(E->ncid, E->vid_lat_rho, &E->lat_rho[0][0])))
        fail("put_var_ failed for lat_rho. Error code = %d\n",E->retval);

    // lon_rho
    if ((E->retval = nc_put_var_double(E->ncid, E->vid_lon_rho, &E->lon_rho[0][0])))
        fail("put_var_ failed. Error code = %d\n",E->retval);

	/*
	// ocean_time
    if ((E->retval = nc_put_var_double(E->ncid, E->vid_ocean_time, &E->interp_time[0])))
        fail("put_var_ failed. Error code = %d\n",E->retval);
		*/

}


void write_coastal_data(e *E) {

	int	ncid;
	int varid;
	int retval;
	size_t attlen = 0;

	int lat_dimid, lon_dimid, time_dimid, time_waves_dimid, time_tides_dimid, dimIds[2], dimIds3d[3], dimIds3d_waves[3], dimIds3d_tides[3] ;
	int lat_varid, lon_varid, wave_time_varid, tide_time_varid, coastline_varid, setup_coast_varid, zeta_coast_varid, tide_coast_varid, added_coast_varid;

	int	Hs_coast_varid, Tp_coast_varid;
	int ocean_time_varid;

	// set the compression parameters
  int shuffle = 1;
  int deflate = 1;
  int deflate_level = 1;

	// create the file
	nc_create(E->fname, NC_NETCDF4|NC_CLOBBER, &ncid);
	// def dimensions
	nc_def_dim(ncid, "time", E->nTimeRoms, &time_dimid);
	//nc_def_dim(ncid, "wave_time", E->nTimeWaves, &time_waves_dimid);
	//nc_def_dim(ncid,"wave_time", E->nTimeWavesSubset, &time_waves_dimid);
	//nc_def_dim(ncid,"tide_time", E->nTimeTideSubset, &time_tides_dimid);
	nc_def_dim(ncid, "lat", E->nLatRho, &lat_dimid);
	nc_def_dim(ncid, "lon", E->nLonRho, &lon_dimid);
	// def vars
	dimIds[0] = lon_dimid;
	dimIds[1] = lat_dimid;

	dimIds3d[0] = time_dimid;
	dimIds3d[1] = lon_dimid;
	dimIds3d[2] = lat_dimid;

	/*
	dimIds3d_waves[0] = time_waves_dimid;
	dimIds3d_waves[1] = lon_dimid;
	dimIds3d_waves[2] = lat_dimid;

	dimIds3d_tides[0] = time_tides_dimid;
	dimIds3d_tides[1] = lon_dimid;
	dimIds3d_tides[2] = lat_dimid;
	*/

	/*
	nc_def_var(ncid, "wave_time", NC_DOUBLE, 1, &dimIds3d_waves[0], &wave_time_varid);
	nc_put_att_text(ncid, wave_time_varid, "long_name", strlen("time"), "time");
	nc_put_att_text(ncid, wave_time_varid, "units", strlen(E->roms_time_units), E->roms_time_units);
	nc_put_att_text(ncid, wave_time_varid, "calendar", strlen("gregorian"), "gregorian");
	*/

	nc_def_var(ncid, "time", NC_DOUBLE, 1, &dimIds3d[0], &ocean_time_varid);
    // define the compression options for this variable
    nc_def_var_deflate(ncid, ocean_time_varid, shuffle, deflate, deflate_level);
		nc_put_att_text(ncid, ocean_time_varid, "long_name", strlen("time"), "time");
		nc_put_att_text(ncid, ocean_time_varid, "units", strlen(E->roms_time_units), E->roms_time_units);
		nc_put_att_text(ncid, ocean_time_varid, "calendar", strlen("gregorian"), "gregorian");
		nc_put_att_text(ncid, ocean_time_varid, "_CoordinateAxisType", strlen("time"), "time");
		nc_put_att_text(ncid, ocean_time_varid, "axis", strlen("T"), "T");


	/*
	nc_def_var(ncid, "tide_time", NC_DOUBLE, 1, &dimIds3d_tides[0], &tide_time_varid);
	nc_put_att_text(ncid, ocean_time_varid, "long_name", strlen("time"), "time");
	nc_put_att_text(ncid, ocean_time_varid, "units", strlen(E->roms_time_units), E->roms_time_units);
	nc_put_att_text(ncid, ocean_time_varid, "calendar", strlen("gregorian"), "gregorian");
	*/

	nc_def_var(ncid, "lat", NC_DOUBLE, 2, dimIds, &lat_varid);
	// define the compression options for this variable
  nc_def_var_deflate(ncid, lat_varid, shuffle, deflate, deflate_level);
	nc_put_att_text(ncid, lat_varid, "long_name", strlen("latitude"), "latitude");
	nc_put_att_text(ncid, lat_varid, "standard_name", strlen("latitude"), "latitude");
	nc_put_att_text(ncid, lat_varid, "units", strlen("degrees north"), "degrees north");
	nc_put_att_text(ncid, lat_varid, "axis", strlen("Y"), "Y");
	//nc_put_att_text(ncid, lat_varid, "_CoordinateAxisType", strlen("Lat"), "Lat");


	nc_def_var(ncid, "lon", NC_DOUBLE, 2, dimIds, &lon_varid);
	// define the compression options for this variable
  nc_def_var_deflate(ncid, lon_varid, shuffle, deflate, deflate_level);
	nc_put_att_text(ncid, lon_varid, "long_name", strlen("Lon"), "longitude");
	nc_put_att_text(ncid, lon_varid, "standard_name", strlen("longitude"), "longitude");
	nc_put_att_text(ncid, lon_varid, "units", strlen("degrees east"), "degrees east");
	nc_put_att_text(ncid, lon_varid, "axis", strlen("X"), "X");
	//nc_put_att_text(ncid, lon_varid, "_CoordinateAxisType", strlen("Lon"), "Lon");



	nc_def_var(ncid, "coastline_mask", NC_DOUBLE, 2, dimIds, &coastline_varid);
	// define the compression options for this variable
  nc_def_var_deflate(ncid, coastline_varid, shuffle, deflate, deflate_level);
	nc_put_att_text(ncid, coastline_varid, "long_name", strlen("land sea mask"), "land sea mask");
	nc_put_att_text(ncid, coastline_varid, "flag_values", strlen("0, 1"), "0, 1");
	nc_put_att_text(ncid, coastline_varid, "flag_meaning", strlen("land water"), "land water");
	nc_put_att_text(ncid, coastline_varid, "coordinates", strlen("lat lon"), "lat lon");



	nc_def_var(ncid, "surge", NC_DOUBLE, 3, dimIds3d, &zeta_coast_varid);
	// define the compression options for this variable
  nc_def_var_deflate(ncid, zeta_coast_varid, shuffle, deflate, deflate_level);
	nc_put_att_text(ncid, zeta_coast_varid, "long_name", strlen("sea level due to atmospheric wind and pressure"), "sea level due to atmospheric wind and pressure");
	nc_put_att_text(ncid, zeta_coast_varid, "standard_name", strlen("sea_surface_height_above_sea_level"), "sea_surface_height_above_sea_level");
	nc_put_att_text(ncid, zeta_coast_varid, "units", strlen("m"), "m");
	nc_put_att_text(ncid, zeta_coast_varid, "positive", strlen("up"), "up");
	nc_put_att_text(ncid, zeta_coast_varid, "coordinates", strlen("lat lon"), "lat lon");


	//nc_def_var(ncid, "setup", NC_DOUBLE, 3, dimIds3d_waves, &setup_coast_varid);
	nc_def_var(ncid, "wave_setup", NC_DOUBLE, 3, dimIds3d, &setup_coast_varid);
	// define the compression options for this variable
  nc_def_var_deflate(ncid, setup_coast_varid, shuffle, deflate, deflate_level);
	nc_put_att_text(ncid, setup_coast_varid, "long_name", strlen("sea level due to wave setup"), "sea level due to wave setup");
	nc_put_att_text(ncid, setup_coast_varid, "standard_name", strlen("sea_surface_height_above_sea_level"), "sea_surface_height_above_sea_level");
	nc_put_att_text(ncid, setup_coast_varid, "units", strlen("m"), "m");
	nc_put_att_text(ncid, setup_coast_varid, "positive", strlen("up"), "up");
	nc_put_att_text(ncid, setup_coast_varid, "method", strlen("simple paramaterisation"), "simple paramaterisation");
	nc_put_att_text(ncid, setup_coast_varid, "coordinates", strlen("lat lon"), "lat lon");



	nc_def_var(ncid, "sea_level", NC_DOUBLE, 3, dimIds3d, &added_coast_varid);
	// define the compression options for this variable
  nc_def_var_deflate(ncid, added_coast_varid, shuffle, deflate, deflate_level);
	nc_put_att_text(ncid, setup_coast_varid, "long_name", strlen("sea level"), "sea level");
	nc_put_att_text(ncid, setup_coast_varid, "standard_name", strlen("sea_surface_height_above_sea_level"), "sea_surface_height_above_sea_level");
	nc_put_att_text(ncid, setup_coast_varid, "units", strlen("m"), "m");
	nc_put_att_text(ncid, setup_coast_varid, "positive", strlen("up"), "up");
	nc_put_att_text(ncid, setup_coast_varid, "method", strlen("linear superposition"), "linear superposition");
	nc_put_att_text(ncid, setup_coast_varid, "ancillary_variables", strlen("level_MSL"), "level_MSL");
	nc_put_att_text(ncid, setup_coast_varid, "coordinates", strlen("lat lon"), "lat lon");


	/*
	nc_def_var(ncid, "Hs", NC_DOUBLE, 3, dimIds3d_waves, &Hs_coast_varid);
	nc_put_att_text(ncid, Hs_coast_varid, "coordinates", strlen("lat_rho lon_rho"), "lat_rho lon_rho");

	nc_def_var(ncid, "Tp", NC_DOUBLE, 3, dimIds3d_waves, &Tp_coast_varid);
	nc_put_att_text(ncid, Tp_coast_varid, "coordinates", strlen("lat_rho lon_rho"), "lat_rho lon_rho");

	*/


	nc_def_var(ncid, "tide", NC_DOUBLE, 3, dimIds3d, &tide_coast_varid);
	// define the compression options for this variable
	nc_def_var_deflate(ncid, tide_coast_varid, shuffle, deflate, deflate_level);
	nc_put_att_text(ncid, tide_coast_varid, "long_name", strlen("sea level due to harmonic tide"), "sea level due to harmonic tide");
	nc_put_att_text(ncid, tide_coast_varid, "standard_name", strlen("sea_surface_height_above_sea_level"), "sea_surface_height_above_sea_level");
	nc_put_att_text(ncid, tide_coast_varid, "units", strlen("m"), "m");
	nc_put_att_text(ncid, tide_coast_varid, "positive", strlen("up"), "up");
	nc_put_att_text(ncid, tide_coast_varid, "coordinates", strlen("lat lon"), "lat lon");


	nc_enddef(ncid);
	// write the data

	nc_put_var_double(ncid, ocean_time_varid, &E->romsTime[0]);
	//nc_put_var_double(ncid, wave_time_varid, &E->wavesTime[0]);
	//nc_put_var_double(ncid, wave_time_varid, &E->waves_interp_time[0]);
	nc_put_var_double(ncid, lat_varid, &E->lat_rho[0][0]);
	nc_put_var_double(ncid, lon_varid, &E->lon_rho[0][0]);
	nc_put_var_double(ncid, coastline_varid, &E->coastline_mask[0][0]);
	nc_put_var_double(ncid, zeta_coast_varid, &E->zeta_coast[0][0][0]);
	//nc_put_var_double(ncid, setup_coast_varid, &E->setup_on_roms[0][0][0]);
	nc_put_var_double(ncid, setup_coast_varid, &E->setup_on_roms_time_interp[0][0][0]);
	//nc_put_var_double(ncid, Hs_coast_varid, &E->Hs_on_roms[0][0][0]);
	//nc_put_var_double(ncid, Tp_coast_varid, &E->Tp_on_roms[0][0][0]);

	nc_put_var_double(ncid, tide_coast_varid, &E->tide_on_roms_time_interp[0][0][0]);

	nc_put_var_double(ncid, added_coast_varid, &E->added[0][0][0]);

	// close the file
	nc_close(ncid);
}



void get_coastal_slope(e *E){

	int	ncid;
	int varid;
	int retval;

    // read in the coastal slope data
    if((retval = nc_open("static_data/slope.nc", NC_NOWRITE, &ncid)))
        fail("failed to open coastal slope data input file: error is %d\n",retval);

	// get the nSlope dimlength
	if((retval = nc_inq_dimid(ncid, "npoints", &varid)))
        fail("failed to get slope npoints dimid: error is %d\n",retval);

    if((retval = nc_inq_dimlen(ncid,varid,&E->nSlopes)))
        fail("failed to get slopes npoints dimlen: error is %d\n",retval);

	// malloc room
	E->slopeLon = malloc(E->nSlopes*sizeof(double));
	E->slopeLat = malloc(E->nSlopes*sizeof(double));
	E->slope = malloc(E->nSlopes*sizeof(double));

	// now get the data
	nc_inq_varid(ncid, "lat", &varid);
    if((retval = nc_get_var_double(ncid, varid, &E->slopeLat[0])))
		fail("failed to read slope lat data: error is %d\n", retval);

	nc_inq_varid(ncid, "lon", &varid);
    if((retval = nc_get_var_double(ncid, varid, &E->slopeLon[0])))
		fail("failed to read slope lon data: error is %d\n", retval);

	nc_inq_varid(ncid, "slope", &varid);
    if((retval = nc_get_var_double(ncid, varid, &E->slope[0])))
		fail("failed to read slope  data: error is %d\n", retval);
}
