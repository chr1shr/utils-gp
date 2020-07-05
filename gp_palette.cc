#include "gp_palette.hh"
#include "gp_pal_defs.hh"

#include <cmath>

#ifdef HAS_PNG
#include <cstdlib>
#include <png.h>

/** Prints an error message due to outputting a PNG image
 * \param[in] err_msg the error message to print. */
void gp_matrix_png_abort(const char* err_msg) {
    fprintf(stderr,"PNG output: %s\n",err_msg);
    exit(1);
}

/** Writes a PNG image using the libpng library.
 * \param[in] (m_,n_) the dimensions of the image.
 * \param[in] rowp a pointer to the PNG bitmap information.
 * \param[in] filename the output file name. */
void gp_matrix_png_write(int m_,int n_,png_bytep *rowp,const char *filename) {
    FILE *fp=fopen(filename,"wb");
    if(fp==NULL) gp_matrix_png_abort("can't open output file");

    png_structp p=png_create_write_struct(PNG_LIBPNG_VER_STRING,NULL,NULL,NULL);
    if(!p) gp_matrix_png_abort("error creating write structure");

    png_infop info=png_create_info_struct(p);
    if(!p) gp_matrix_png_abort("error creating info structure");

    if(setjmp(png_jmpbuf(p))) gp_matrix_png_abort("setjmp error");
    png_init_io(p,fp);
    png_set_IHDR(p,info,m_,n_,8,PNG_COLOR_TYPE_RGB,PNG_INTERLACE_NONE,
                 PNG_COMPRESSION_TYPE_DEFAULT,PNG_FILTER_TYPE_DEFAULT);

    png_write_info(p,info);
    png_write_image(p,rowp);
    png_write_end(p,NULL);
    fclose(fp);
    png_destroy_write_struct(&p,&info);
}
#endif

/** Converts a scalar value into an (R,G,B) color. This internal routine
 * assumes that the scalar value is already scaled and clipped into the [0,1]
 * range.
 * \param[in] fv the scaled and clipped scalar.
 * \param[in] (re,gr,bl) the components of the color. */
void gp_palette::col_internal(float fv,float &re,float &gr,float &bl) {
    const float *cp,*cp2;
    if(equal) {

        // For an equally-spaced palette, use an integer conversion to
        // determine which interval the value is in
        int k=static_cast<int>(fv*(n-1));
        if(k<0) k=0;else if(k>(n-1)) k=n-1;
        fv=fv*(n-1)-k;
        cp=c+3*k;cp2=cp+3;
    } else {

        // For an unequally-spaced palette, step through the spacings to find
        // which interval the value is in
        cp=c;
        while(cp<c+4*(n-1)&&fv>cp[7]) cp+=4;
        fv=(fv-cp[3])/(cp[7]-cp[3]);
        cp2=cp+4;
    }

    // Calculate the (R,G,B) components of the color using linear interpolation
    re=*cp*(1-fv)+*cp2*fv;
    gr=cp[1]*(1-fv)+cp2[1]*fv;
    bl=cp[2]*(1-fv)+cp2[2]*fv;
}

/** Writes the palette in a format that can be read as a Gnuplot command.
 * \param[in] fp the file handle to write to.
 * \param[in] gamma the gamma value to apply to the palette. */
void gp_palette::output_gnuplot(FILE *fp,double gamma) {

    // Write the command and the first color
    fprintf(fp,"set palette defined (0. \"#%06x\"",hex_col(c));

    // Write the colors, using separate approach equal and unequal palettes
    if(equal) {
        double fac=1./(n-1);
        const float *cp=c+3;
        for(int k=1;k<n;k++,cp+=3) fprintf(fp,",%g \"#%06x\"",pow(k*fac,gamma),hex_col(cp));
    } else for(const float *cp=c+4;cp<c+4*n;cp+=4)
        fprintf(fp,",%g \"#%06x\"",pow(cp[3],gamma),hex_col(cp));
    fputs(")\n",fp);
}

/** Writes the palette as a simple text file, which can be used for graphing
 * and analysis.
 * \param[in] fp the file handle to write to.
 * \param[in] gamma the gamma value to apply to the palette. */
void gp_palette::output_text(FILE *fp,double gamma) {
    if(equal) {
        double fac=1./(n-1);
        const float *cp=c;
        for(int k=0;k<n;k++,cp+=3) fprintf(fp,"%g %g %g %g\n",pow(k*fac,gamma),*cp,cp[1],cp[2]);
    } else for(const float *cp=c;cp<c+4*n;cp+=4)
        fprintf(fp,"%g %g %g %g\n",cp[3],*cp,cp[1],cp[2]);
}

/** Outputs a PNG color bar key of the palette used in the bitmap image output.
 * \param[in] filename the name of the file to write to.
 * \param[in] (wid,hei) the dimensions of the image.
 * \param[in] gamma the gamma value to apply to the color bar. */
void gp_palette::colorbar_horizontal(const char* filename,int wid,int hei,double gamma) {
#ifdef HAS_PNG

    // Create PNG bitmap array image. Since all rows are the same, set all the
    // pointers to the zeroth entry.
    png_bytep *rowp=new png_bytep[hei];
    png_byte* rp=(*rowp=new png_byte[3*wid]);
    for(int j=1;j<hei;j++) rowp[j]=rp;

    // Assemble the color bar
    float iwid=1.0/float(wid);
    for(int i=0;i<wid;i++) col(pow((0.5+i)*iwid,gamma),rp);

    // Output the PNG image and free memory
    gp_matrix_png_write(wid,hei,rowp,filename);
    delete [] *rowp;
    delete [] rowp;
#else

    // Give a status message in the case where the code was compiled without
    // the PNG library
    puts("No PNG support");
#endif
}

/** Outputs a PNG color bar key of the palette used in the bitmap image output.
 * \param[in] filename the name of the file to write to.
 * \param[in] (wid,hei) the dimensions of the image.
 * \param[in] gamma the gamma value to apply to the color bar. */
void gp_palette::colorbar_vertical(const char* filename,int wid,int hei,double gamma) {
#ifdef HAS_PNG

    // Create PNG bitmap array image. Since all rows are the same, set all the
    // pointers to the zeroth entry.
    png_bytep *rowp=new png_bytep[hei];

    // Assemble the color bar
    float ihei=1.0/float(hei);
    for(int j=0;j<hei;j++) {
        png_byte* rp=(rowp[j]=new png_byte[3*wid]);
        col(pow((hei-0.5-j)*ihei,gamma),rp);
        for(;rp<rowp[j]+3*wid;rp+=3) {
            *rp=*(rowp[j]);
            rp[1]=rowp[j][1];
            rp[2]=rowp[j][2];
        }
    }

    // Output the PNG image and free memory
    gp_matrix_png_write(wid,hei,rowp,filename);
    for(int j=0;j<hei;j++) delete rowp[j];
    delete [] rowp;
#else

    // Give a status message in the case where the code was compiled without
    // the PNG library
    puts("No PNG support");
#endif
}
