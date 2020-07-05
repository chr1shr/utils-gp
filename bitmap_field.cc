#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>

#include "gp_matrix.hh"
#include "gp_pal_defs.hh"

void syntax_message() {
    fputs("Usage: ./bitmap_field [-g <gamma>] [-h] [-s] [-p <palette>] [-u <upsample>]\n"
          "                      [-c <ilo> <ihi> <jlo> <jhi>]\n"
          "                      <field> <output_png> [<zlo> <zhi>]\n",stderr);
    exit(1);
}

int main(int argc,char **argv) {

    // Check for sufficient command-line arguments
    if(argc<3) syntax_message();

    // Check for the -p option and if present read the exponent
    float gamma=1.;int ups=1,i=1,p_ind=0,ilo=0,ihi=0,jlo=0,jhi=0,h=0;
    bool stagger=false,power=false,crop=false;
    while(i+2<argc) {
        if(strcmp(argv[i],"-g")==0) {
            power=true;gamma=atof(argv[i+1]);i+=2;
        } else if(strcmp(argv[i],"-p")==0) {
            p_ind=gp_pal_index(argv[i+1]);
            if(p_ind==-1) {
                fputs("Unknown palette type\n",stderr);
                return 1;
            }
            i+=2;
        } else if(strcmp(argv[i],"-u")==0) {
            ups=atof(argv[i+1]);
            if(ups<=0||ups>32) {
                fputs("Error in upsampling factor\n",stderr);
                return 1;
            }
            i+=2;
        } else if(strcmp(argv[i],"-s")==0) {
            stagger=true;
            i++;
        } else if(strcmp(argv[i],"-h")==0) {
            h++;i++;
        } else if(strcmp(argv[i],"-c")==0) {
            crop=true;
            if(i+6>argc) {
                fputs("Not enough arguments for cropping option\n",stderr);
                return 1;
            }
            ilo=atoi(argv[++i]);ihi=atoi(argv[++i]);
            jlo=atoi(argv[++i]);jhi=atoi(argv[++i]);
            i++;
            if(ilo<0) {fprintf(stderr,"ilo mapped from %d to 0\n",ilo);ilo=0;}
            if(jlo<0) {fprintf(stderr,"jlo mapped from %d to 0\n",jlo);jlo=0;}
        } else break;
    }

    // Check for the correct number of command-line arguments
    if(argc!=i+4&&argc!=i+2) syntax_message();

    // Read in the field and apply downsampling
    gp_matrix gp(argv[i]);
    while(h>0) {gp.half_downsample();h--;}

    // Apply cropping and the power transform if needed
    if(crop) {
        if(ihi>gp.m) {fprintf(stderr,"ihi mapped from %d to %d\n",ihi,gp.m);ihi=gp.m;}
        if(jhi>gp.n) {fprintf(stderr,"jhi mapped from %d to %d\n",jhi,gp.n);jlo=gp.n;}
        gp.crop(ilo,ihi,jlo,jhi);
    }
    if(power) gp.power(gamma);

    // Calculate the z scale, or use the specified values
    float zlo,zhi;
    if(argc==i+4) zlo=atof(argv[i+2]),zhi=atof(argv[i+3]);
    else {
        gp.field_range(zlo,zhi);
        printf("%g %g\n",zlo,zhi);
    }

    // Output the bitmap
    if(ups>1) gp.output_bitmap_upsample(argv[i+1],ups,stagger,zlo,zhi,p_ind);
    else gp.output_bitmap(argv[i+1],zlo,zhi,p_ind);
}
