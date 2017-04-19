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

void write_json(e *E){

    FILE *o;
    char fname[80];
    int t,i,j;
    int count = 0;

    sprintf(fname,"storm_coast.json");
    o = fopen(fname,"w");


    fprintf(o,"{\"type\": \"FeatureCollection\",\"properties\":");
    fprintf(o,"{\"basedatetime\":\"%.0f\"},", E->romsTime[0]);
    fprintf(o,"\"features\":");
    fprintf(o,"[{\"geometry\":");

    count = 0;
    for(i=0;i<E->nLonRho;i++){
      for(j=0;j<E->nLatRho;j++){
        if(E->setup_on_roms_time_interp[0][i][j] != NC_FILL_DOUBLE){
          fprintf(o,"{\"type\":");
          fprintf(o,"\"Point\",");
          fprintf(o,"\"coordinates\":[\"%0.6f\",\"%0.6f\"]},", E->lon_rho[i][j],E->lat_rho[i][j]);
          fprintf(o,"\"type\": \"Feature\",\"properties\":");
          fprintf(o,"{\"station\": \"%010d\"", count);
          fprintf(o,"}},");
          count++;
        }
      }
    }
    fprintf(o,"]}\n");

    fclose(o);

    // write point location info
/*

    {"type": "FeatureCollection", 
  "properties": 
  {"basedatetime": "201703101200"}, 
  "features": 
    [{"geometry": 
      {"type": 
        "Point", 
        "coordinates": ["113.397499", "-24.220100"]}, 
        "type": "Feature", 
        "properties": 
          {"bias_set": "True", 
          "abs min forecast": "0.54", 
          "alert": "False", 
          "forecast exceedance of HAT": "-0.37", 
          "station": "006108", 
          "abs max forecast": "1.48", 
          "station_name": "CAPE CUVIER WHARF"
          }
      }, 
      {"geometry": 
           {"type": "Point", 
            "coordinates": ["152.400604", "-24.761101"]
           }, 
           "type": "Feature", 
           "properties": 
            {"bias_set": "True", 
            "abs min forecast": "0.58", 
            "alert": "False", 
            "forecast exceedance of HAT": "-0.32", 
            "station": "539076", 
            "abs max forecast": "3.33", 
            "station_name": "BURNETT HEADS TIDE TM"
            }
        } 
      ]}

*/



    // {"type": "FeatureCollection", "properties": {"time_regular_inc": 3600, "basedatetime": "201603180000", "features": [{"geometry": 
    //
    // {"type": "Point", 
    //     "coordinates": ["122.218300", "-18.000799"]}, 
    //         "type": "Feature", 
    //                 "properties": 
    //                             {"eta_fc": 
    //                                             [["1458261000", "0.160415"],



    // write time series location info for this location
    // write header 
    //fprintf(out,"{\"type\": \"FeatureCollection\", \"properties\": {\"time_regular_inc": %d, \"basedatetime\": \"201603180000\", \"features\": [{\"geometry\":"); // need to add time interval, basedatetime




}
