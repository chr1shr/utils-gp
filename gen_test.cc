#include <cstdio>
#include <cstdlib>
#include <cmath>

#include "gp_matrix.hh"

int main(int argc,char **argv) {

    // Check for the correct number of command-line arguments
    if(argc!=2) {
        fputs("Usage: ./gen_test <grid_size>\n",stderr);
        return 1;
    }

    // Read in the grid resolution and check it is in a reasonable range
    int m=atoi(argv[1]);
    if(m<=0||m>65536) {
        fputs("Grid size out of range\n",stderr);
        return 1;
    }

    // Create the Gnuplot matrix class with the given grid resolution, and set
    // the coordinate ranges
    gp_matrix fld(m,m);
    double s=0.5/m;
    fld.set_coordinate_ranges(s,1-s,s,1-s);

    // Fill the grid with the waves pattern, and output it
    float *xp,*yp,*fp=fld.f,x,y;
    for(yp=fld.y;yp<fld.y+m;yp++) {
        y=*yp;
        for(xp=fld.x;xp<fld.x+m;xp++,fp++)
            x=*xp,
            *fp=sin((y+exp(2*x)+0.4*sin(20*y+exp(2.5*x)))*5)+0.3*cos(40*x*y);
    }
    fld.output("waves.fe");

    // Fill the grid with the superellipse pattern, and output it
    for(fp=fld.f,yp=fld.y;yp<fld.y+m;yp++) {
        y=*yp-0.5;y*=y;y*=y;
        for(xp=fld.x;xp<fld.x+m;xp++,fp++)
            x=*xp-0.5,x*=x,x*=x,
            *fp=sqrt(sqrt(x+y))-0.4;
    }
    fld.output("sellipse.fe");

}
