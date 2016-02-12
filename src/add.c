#include "grid.h"

void add(e *E){

    int t,i,j;

    E->added = malloc3d_double(E->nTimeRoms, E->nLonRho, E->nLatRho);

    // initialize this array
    for(t=0;t<E->nTimeRoms;t++){
        for(i=0;i<E->nLonRho;i++){
		    for(j=0;j<E->nLatRho;j++){
                E->added[t][i][j] = NC_FILL_DOUBLE;
            }
        }
    }

    // now add them
    for(t=0;t<E->nTimeRoms;t++){
        for(i=0;i<E->nLonRho;i++){
		    for(j=0;j<E->nLatRho;j++){
                if(E->setup_on_roms_time_interp[t][i][j] != NC_FILL_DOUBLE){
                    E->added[t][i][j] = E->setup_on_roms_time_interp[t][i][j] +
                                        E->tide_on_roms_time_interp[t][i][j] +
                                        E->zeta_coast[t][i][j];
                }
            }
        }
    }

}

void write_time_series(e *E){

    FILE    *out;
    char    fname[80];
    int     t,i,j;

    for(i=0;i<E->nLonRho;i++){
	    for(j=0;j<E->nLatRho;j++){
            if(E->setup_on_roms_time_interp[0][i][j] != NC_FILL_DOUBLE){
                sprintf(fname,"trash/%.4f_%.4f.txt",E->lon_rho[i][j], E->lat_rho[i][j]);
                out = fopen(fname,"w");
                fprintf(out,"# zeta\ttide\tsetup\ttotal\n");
                for(t=0;t<E->nTimeRoms;t++){
                    fprintf(out,"%.4f %.4f %.4f %.4f\n", E->zeta_coast[t][i][j], E->tide_on_roms_time_interp[t][i][j],
                    E->setup_on_roms_time_interp[t][i][j],
                    E->added[t][i][j]);
                }
                fclose(out);
            }
        }
    }

}
