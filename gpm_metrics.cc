#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "gp_matrix.hh"

int main(int argc,char **argv) {

    // Check for the correct number of command-line arguments
    if(argc<2) {
        fputs("Usage: ./gpm_metrics [-d] [-f] [-s] <input_field>\n\n"
              "Reads in a gnuplot matrix file and calculates various metrics. If no options\n"
              "are specified, then all metrics are printed in a human-readable format.\n"
              "Otherwise, specific metrics are printed as numbers for automatic parsing.\n\n"
              "Options: -d the dimensions of the matrix\n"
              "         -e the minimum and maximum field values\n"
              "         -a the average field value\n"
              "         -s the estimated grid spacings\n",stderr);
        return 1;
    }

    // Parse the command-line input for output flags
    unsigned int oflags=0;
    for(int k=1;k<argc-1;k++) {
        if(strcmp(argv[k],"-d")==0) oflags|=1;
        else if(strcmp(argv[k],"-e")==0) oflags|=2;
        else if(strcmp(argv[k],"-a")==0) oflags|=4;
        else if(strcmp(argv[k],"-s")==0) oflags|=8;
        else {
            fputs("Error reading command-line arguments\n",stderr);
            return 1;
        }
    }

    // Read in the Gnuplot file. If it doesn't exist, there will be an error.
    gp_matrix fld(argv[argc-1]);
    int &m=fld.m,&n=fld.n;
    if(oflags==0) {

        // If no options were provided, then print the human-readable
        // information
        float zlo,zhi;
        fld.field_range(zlo,zhi);
        printf("%s: %d by %d points, (min,avg,max)=(%g,%g,%g) (dx,dy)=(%g,%g)\n",
            argv[argc-1],m,n,zlo,fld.average(),zhi,
            (fld.x[m-1]-*fld.x)/(m-1),(fld.y[n-1]-*fld.y)/(n-1));
    } else {

        // If some options were provided, then print the corresponding
        // information as plain text numbers
        if(oflags&1) printf("%d %d%c",m,n,oflags&14?' ':'\n');
        if(oflags&2) {
            float zlo,zhi;
            fld.field_range(zlo,zhi);
            printf("%g %g%c",zlo,zhi,oflags&12?' ':'\n');
        }
        if(oflags&4) printf("%g%c",fld.average(),oflags&8?' ':'\n');
        if(oflags&8) printf("%g %g\n",(fld.x[m-1]-*fld.x)/(m-1),(fld.y[n-1]-*fld.y)/(n-1));
    }
}
