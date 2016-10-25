#include "grid.h"

void get_cli_args(e *E, int argc, char *argv[]){

  int             c;
  const char    * short_opt = "r:t:w:o:l:h";
  struct option   long_opt[] =
  {
     {"help",          no_argument,       NULL, 'h'},
     {"roms",          required_argument, NULL, 'r'},
     {"tides",          required_argument, NULL, 't'},
     {"waves",          required_argument, NULL, 'w'},
     {"output",          required_argument, NULL, 'o'},
     {"levels",          required_argument, NULL, 'l'},
     {NULL,            0,                 NULL, 0  }
  };

  E->haveRoms = 0;
 	E->haveTides = 0;
  E->haveWaves = 0;
  E->haveLevels = 0;
  E->haveOutput = 0;

  while((c = getopt_long(argc, argv, short_opt, long_opt, NULL)) != -1)
  {
     switch(c)
     {
        case -1:       /* no more arguments */
        case 0:        /* long options toggles */
         break;

        case 'r':
          E->roms_input = malloc((strlen(optarg)+1)*sizeof(char));
          strncpy(E->roms_input, optarg, strlen(optarg));
          // fix the string
          E->roms_input[strlen(optarg)] = '\x0';
          printf("roms_input file: %s\n", E->roms_input);
          E->haveRoms = 1;
          break;

        case 't':
          E->tide_input = malloc((strlen(optarg)+1)*sizeof(char));
          strncpy(E->tide_input, optarg, strlen(optarg));
          // fix the string
          E->tide_input[strlen(optarg)] = '\x0';
          printf("tides_file = %s\n", E->tide_input);
          E->haveTides = 1;
          break;

        case 'w':
          E->wave_input = malloc((strlen(optarg)+1)*sizeof(char));
          strncpy(E->wave_input, optarg, strlen(optarg));
          // fix the string
          E->wave_input[strlen(optarg)] = '\x0';
          printf("waves_file = %s\n", E->wave_input);
          E->haveWaves = 1;
          break;

        case 'o':
          E->fname = malloc((strlen(optarg)+1)*sizeof(char));
          strncpy(E->fname, optarg, strlen(optarg));
          // fix the string
          E->fname[strlen(optarg)] = '\x0';
          printf("output_file = %s\n", E->fname);
          E->haveOutput = 1;
          break;

        case 'l':
          E->levels_input = malloc((strlen(optarg)+1)*sizeof(char));
          strncpy(E->levels_input, optarg, strlen(optarg));
          // fix the string
          E->levels_input[strlen(optarg)] = '\x0';
          printf("levels_file = %s\n", E->levels_input);
          E->haveLevels= 1;
          break;

        case 'h':
          print_usage();
          exit(0);

        default:
          print_usage();
          exit(0);
     };
  };

  // check presence
  if(  E->haveRoms == FALSE  ){
    fprintf(stderr,"no reference ROMS grid file specified:\nI need one of these to run.\n");
    exit(1);
  }
  if(  E->haveOutput == FALSE ){
    fprintf(stderr,"no output file specified:\nI need one of these to run.\n");
    exit(1);
  }
  if(  E->haveWaves == FALSE ){
    fprintf(stderr,"no wave input file specified:\nI need one of these to run.\n");
    exit(1);
  }
  if(  E->haveTides == FALSE ){
    fprintf(stderr,"no tides input file specified:\nI need one of these to run.\n");
    exit(1);
  }
  //if(  E->haveLevels == FALSE ){
  //  fprintf(stderr,"no reference level input file specified:\nI need one of these to run.\n");
  //  exit(1);
  //}
}


void print_usage(){
  printf("Usage: \n");
  printf("  -h, --help                print this help and exit\n");
  printf("  -r, --roms [filename]            roms input netcdf file\n");
  printf("  -t, --tides [filename]           tides input netcdf file\n");
  printf("  -w, --waves [filename]           waves input netcdf file\n");
  printf("  -l, --levels [filename]          reference levels input netcdf file\n");
  printf("  -o, --output [filename]          output netcdf file\n");
  printf("\n");
}
