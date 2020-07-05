#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>

#include "gp_pal_defs.hh"
#include "gp_palette.hh"

int main(int argc,char **argv) {

    // Check for the correct number of command-line arguments
    if(argc<6||argc>7) {
        fputs("Usage: ./colorbar {h|v} <palette> <width> <height> <output_png> [gamma]\n",stderr);
        return 1;
    }

    // Check the dimensions are sensible
    int wid=atoi(argv[3]),hei=atoi(argv[4]);
    if(wid<=0||wid>16777216||hei<=0||hei>16777216) {
        fputs("Image dimensions out of bounds\n",stderr);
        return 1;
    }

    // Find the palette type
    int p_ind=gp_pal_index(argv[2]);
    if(p_ind==-1) {
        fputs("Unknown palette type; available palettes are:\n\n",stderr);
        gp_list_palettes(stderr);
        return 1;
    }

    // Read the gamma value if specified
    double gamma=1.;
    if(argc==7) {
        gamma=atof(argv[6]);
        if(gamma<=0||gamma>1e6) {
            fputs("Gamma value is outside a typical range\n",stderr);
            return 1;
        }
    }

    // Create and output the color bar
    if(strcmp(argv[1],"h")==0) gp_pals[p_ind].colorbar_horizontal(argv[5],wid,hei,gamma);
    else if(strcmp(argv[1],"v")==0) gp_pals[p_ind].colorbar_vertical(argv[5],wid,hei,gamma);
    else {
        fputs("Mode should be either 'h' (horizontal) or 'v' (vertical)\n",stderr);
        return 1;
    }
}
