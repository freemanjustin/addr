#include "grid.h"

void process_reference_levels(e *E){

  int	ncid;
  int varid;
  int retval;

  // kd tree variables
  double  pos[2];
  int     i,j;
	int count;
	grid_index		*idx;

  // read in the reference levels netcdf data
  // open the file
  if((retval = nc_open(E->levels_input, NC_NOWRITE, &ncid)))
    fail("failed to open reference levels input file: error is %d\n",retval);

	// get the lat dimension sizes
  if((retval = nc_inq_dimid(ncid, "xi_rho", &varid)))
    fail("failed to get roms lat dimid: error is %d\n",retval);

  if((retval = nc_inq_dimlen(ncid,varid,&E->levelXiRho)))
    fail("failed to get reference levels lat dimid: error is %d\n",retval);

  //printf("levels xi_rho = %zu\n", E->levelXiRho);

  // get the lon dimension sizes
  if((retval = nc_inq_dimid(ncid, "eta_rho", &varid)))
    fail("failed to get reference levels lon_rho dimid: error is %d\n",retval);

  if((retval = nc_inq_dimlen(ncid,varid,&E->levelEtaRho)))
    fail("failed to read roms lon_rho data: error is %d\n",retval);

  //printf("levels eta_rho = %zu\n", E->levelEtaRho);

  // malloc room for the arrays

  E->level_lon_rho = malloc2d_double(E->levelEtaRho, E->levelXiRho);
  E->level_lat_rho = malloc2d_double(E->levelEtaRho, E->levelXiRho);
  E->level_AHD = malloc2d_double(E->levelEtaRho, E->levelXiRho);
  E->level_GDA94 = malloc2d_double(E->levelEtaRho, E->levelXiRho);
  E->level_HAT = malloc2d_double(E->levelEtaRho, E->levelXiRho);
  E->level_LAT = malloc2d_double(E->levelEtaRho, E->levelXiRho);
  E->level_MSL = malloc2d_double(E->levelEtaRho, E->levelXiRho);

  nc_inq_varid(ncid, "lat_rho", &varid);
  if((retval = nc_get_var_double(ncid, varid, &E->level_lat_rho[0][0])))
  fail("failed to read reference level lat_rho data: error is %d\n", retval);

  nc_inq_varid(ncid, "lon_rho", &varid);
  if((retval = nc_get_var_double(ncid, varid, &E->level_lon_rho[0][0])))
  fail("failed to read reference level lon_rho data: error is %d\n", retval);

  nc_inq_varid(ncid, "level_AHD", &varid);
  if((retval = nc_get_var_double(ncid, varid, &E->level_AHD[0][0])))
  fail("failed to read reference level level_AHD data: error is %d\n", retval);

  nc_inq_varid(ncid, "level_GDA94", &varid);
  if((retval = nc_get_var_double(ncid, varid, &E->level_GDA94[0][0])))
  fail("failed to read reference level level_GDA94 data: error is %d\n", retval);

  nc_inq_varid(ncid, "level_HAT", &varid);
  if((retval = nc_get_var_double(ncid, varid, &E->level_HAT[0][0])))
  fail("failed to read reference level level_HAT data: error is %d\n", retval);

  nc_inq_varid(ncid, "level_LAT", &varid);
  if((retval = nc_get_var_double(ncid, varid, &E->level_LAT[0][0])))
  fail("failed to read reference level level_LAT data: error is %d\n", retval);

  nc_inq_varid(ncid, "level_MSL", &varid);
  if((retval = nc_get_var_double(ncid, varid, &E->level_MSL[0][0])))
  fail("failed to read reference level level_MSL data: error is %d\n", retval);

  nc_close(ncid);

  /*
  float level_AHD(eta_rho, xi_rho) ;
		level_AHD:_FillValue = -9999.f ;
		level_AHD:long_name = "tidal reference level AHD" ;
		level_AHD:standard_name = "water_surface_height_above_reference_datum" ;
		level_AHD:units = "m" ;
		level_AHD:ancillary_variables = "level_MSL" ;
		level_AHD:method = "manual application of VDT tool to model gridpoints. http://www.crcsi.com.au/research/commissioned-research/auscoast-vdt/" ;
	float level_GDA94(eta_rho, xi_rho) ;
		level_GDA94:_FillValue = -9999.f ;
		level_GDA94:long_name = "tidal reference level GDA94" ;
		level_GDA94:standard_name = "water_surface_height_above_reference_datum" ;
		level_GDA94:units = "m" ;
		level_GDA94:ancillary_variables = "level_MSL" ;
		level_GDA94:method = "manual application of VDT tool to model gridpoints. http://www.crcsi.com.au/research/commissioned-research/auscoast-vdt/" ;
	float level_HAT(eta_rho, xi_rho) ;
		level_HAT:_FillValue = -9999.f ;
		level_HAT:long_name = "tidal reference level HAT" ;
		level_HAT:standard_name = "water_surface_height_above_reference_datum" ;
		level_HAT:units = "m" ;
		level_HAT:ancillary_variables = "level_MSL" ;
		level_HAT:method = "manual application of VDT tool to model gridpoints. http://www.crcsi.com.au/research/commissioned-research/auscoast-vdt/" ;
	float level_LAT(eta_rho, xi_rho) ;
		level_LAT:_FillValue = -9999.f ;
		level_LAT:long_name = "tidal reference level LAT" ;
		level_LAT:standard_name = "water_surface_height_above_reference_datum" ;
		level_LAT:units = "m" ;
		level_LAT:ancillary_variables = "level_MSL" ;
		level_LAT:method = "manual application of VDT tool to model gridpoints. http://www.crcsi.com.au/research/commissioned-research/auscoast-vdt/" ;
	float level_MSL(eta_rho, xi_rho) ;
		level_MSL:_FillValue = -9999.f ;
		level_MSL:long_name = "tidal reference level MSL" ;
		level_MSL:standard_name = "water_surface_height_above_reference_datum" ;
		level_MSL:units = "m" ;
		level_MSL:ancillary_variables = "level_MSL" ;
		level_MSL:method = "manual application of VDT tool to model gridpoints. http://www.crcsi.com.au/research/commissioned-research/auscoast-vdt/" ;

    */

  // we should already have the target roms grid in memory
  // including the coasline mask


  // construct the kdtree for grid searching
  // this is the reference level data
  E->kd = kd_create(2);
  E->roms_index = malloc(E->levelEtaRho*E->levelXiRho*sizeof(grid_index));

  // malloc rooms for the output data on the roms grid
  E->level_AHD_onRoms = malloc2d_double(E->nLonRho, E->nLatRho);
  E->level_GDA94_onRoms = malloc2d_double(E->nLonRho, E->nLatRho);
  E->level_HAT_onRoms = malloc2d_double(E->nLonRho, E->nLatRho);
  E->level_LAT_onRoms = malloc2d_double(E->nLonRho, E->nLatRho);
  E->level_MSL_onRoms = malloc2d_double(E->nLonRho, E->nLatRho);

  count = 0;
  for(i=0;i<E->levelEtaRho;i++){
    for(j=0;j<E->levelXiRho;j++){
      E->roms_index[count].i = i;
      E->roms_index[count].j = j;
      pos[0] = E->level_lon_rho[i][j];
      pos[1] = E->level_lat_rho[i][j];
      kd_insert(E->kd, pos, &E->roms_index[count]);
      count++;
    }
  }

  //printf("done kd_create and kd_inset\n");


  // for each coastal grid cell in the target grid
  // search the kd-tree for the corresponding level data

  count = 0;
	for(i=0;i<E->nLonRho;i++){
		for(j=0;j<E->nLatRho;j++){
      if(E->coastline_mask[i][j] == 1) { // if this point is a coastal point
  		 	pos[0] = E->lon_rho[i][j];
  			pos[1] = E->lat_rho[i][j];
  			E->set = kd_nearest_range(E->kd, pos, 0.02);	// make the .05 dynamic
  			//printf("##point %d ## search around grid point %f, %f returned %d grid items\n", count,pos[0], pos[1], kd_res_size(E->set));
        if(kd_res_size(E->set) == 0) printf("########### missed point!! increase search range ############\n");
        while( !kd_res_end( E->set ) ) {
  		    // get the data and position of the current result item
  		    idx = (grid_index*)kd_res_item( E->set, pos );

  		    // compute the distance of the current result from the pt
  		    //dist = sqrt( dist_sq( pt, pos, 2 ) );

  		    // print out the retrieved data
  		    //printf( "\tpos at (%.3f, %.3f) idx = %d %d\n", pos[0], pos[1], idx->i, idx->j );
          //printf( "\t\tdata at (%.3f, %.3f) [idx = %d %d] = %f\n", pos[0], pos[1], idx->i, idx->j, E->level_AHD[idx->i][idx->j] );

          // this is where we construct the output data for writing
          // by copying the reference level data to the
          // subset arrays

          E->level_AHD_onRoms[i][j] = E->level_AHD[idx->i][idx->j];
          E->level_GDA94_onRoms[i][j] = E->level_GDA94[idx->i][idx->j];
          E->level_HAT_onRoms[i][j] = E->level_HAT[idx->i][idx->j];
          E->level_LAT_onRoms[i][j] = E->level_LAT[idx->i][idx->j];
          E->level_MSL_onRoms[i][j] = E->level_MSL[idx->i][idx->j];

  		    // go to the next entry
  		    kd_res_next( E->set );
  		  }
  			count++;
      }
      else{ // make it a fill value
        E->level_AHD_onRoms[i][j] = NC_FILL_DOUBLE;
        E->level_GDA94_onRoms[i][j] = NC_FILL_DOUBLE;
        E->level_HAT_onRoms[i][j] = NC_FILL_DOUBLE;
        E->level_LAT_onRoms[i][j] = NC_FILL_DOUBLE;
        E->level_MSL_onRoms[i][j] = NC_FILL_DOUBLE;
      }
 	 		//printf("count = %d (outof %d total)\n", count, E->number_of_points_inside);
		}
	}
}
