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
    
    
	// malloc the work struct
	E = malloc(sizeof(e));
	if(E==NULL) fail("Malloc failed\n");
	
	// parse command line arguments
	if(argc < 3){
		fail("need an input xml file and an output filename\n");
	}
	else{
		get_command_line_arg_as_string(&E->input_xml, argv[1]);
		get_command_line_arg_as_string(&E->fname, argv[2]);
	}
	

    // read input XML
    get_params(E);
    
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
    
