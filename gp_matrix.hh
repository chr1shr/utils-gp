#ifndef GP_MATRIX_HH
#define GP_MATRIX_HH

#include <cstdio>
#include <cstdlib>
#include <cmath>

FILE* gp_safe_fopen(const char *filename,const char *mode);

class gp_matrix {
    public:
        /** The horizontal size of the grid. */
        int m;
        /** The vertical size of the grid. */
        int n;
        /** The total number of grid points. */
        int mn;
        /** The x coordinates of the grid. */
        float *x;
        /** The y coordinates of the grid. */
        float *y;
        /** The field values at the grid points. */
        float *f;
        /** The indicator field for constructing contours. */
        int *s;
        /** The horizontal grid spacing, assuming that the grid x coordinates
         * are equally spaced. */
        float dx;
        /** The vertical grid spacing, assuming that the grid y coordinates are
         * equally spaced. */
        float dy;
        /** The reciprocal of the horizontal grid spacing. */
        float xsp;
        /** The reciprocal of the vertical grid spacing. */
        float ysp;
        gp_matrix(const char* filename);
        gp_matrix(int m_,int n_);
        gp_matrix(gp_matrix& gp);
        gp_matrix(gp_matrix &gp,float xmin,float xmax,float ymin,float ymax);
        ~gp_matrix();
        void read(const char* filename);
        void set_coordinate_ranges(double ax,double bx,double ay,double by);
        /** Returns a bilinear interpolation of the field at a given location.
         * \param[in] (xs,ys) the position at which to carry out the
         *                    interpolation.
         * \return The interpolated value. */
        inline float field_lp(float xs,float ys) {
            return field_lp_int((xs-*x)*xsp,(ys-*y)*ysp);
        }
        float field_lp_int(float xi,float yi);
        void x_transect(float ys,float *xt);
        void y_transect(float xs,float *yt);
        void y_extrema(float xs,float &fmin,float &fmax);
        float circle_average(double x,double y,double r);
        float circle_average(double x,double y,double r,double &mi,double &mx);
        void crop(int ilo,int ihi,int jlo,int jhi);
        void power(float gamma);
        float average();
        void zero_project();
        void half_downsample();
        /** Outputs the data in the Gnuplot matrix binary format.
         * \param[in] filename the name of the file to write to. */
        inline void output(const char* filename) {
            FILE *fp=gp_safe_fopen(filename,"wb");
            output(fp);
            fclose(fp);
        }
        void output(FILE *fp=stdout);
        /** Outputs a double-resolution file where each grid point is converted into
         * a rectangular step.
         * \param[in] filename the name of the file to write to. */
        inline void output_stepped(const char *filename) {
            FILE *fp=gp_safe_fopen(filename,"wb");
            output_stepped(fp);
            fclose(fp);
        }
        void output_stepped(FILE *fp=stdout);
        /** Outputs the contours of the field, tracing out each contour as a
         * contiguous line for plotting in Gnuplot.
         * \param[in] hm a pointer to an array of contour values to consider.
         * \param[in] hc the number of elements in the array.
         * \param[in] filename the name of the file to write to. */
        inline void output_contours(float *hm,int hc,const char *filename) {
            FILE *fp=gp_safe_fopen(filename,"w");
            output_contours(hm,hc,fp);
            fclose(fp);
        }
        void output_contours(float *hm,int hc,FILE *fp=stdout);
        /** Outputs the contours of the field, tracing out each contour as a
         * contiguous line for plotting in Gnuplot, and removing parts of the
         * contours in regions where a given level set function is positive.
         * \param[in] hm a pointer to an array of contour values to consider.
         * \param[in] hc the number of elements in the array.
         * \param[in] phi a level set function in a gp_matrix class.
         * \param[in] filename the name of the file to write to. */
        inline void output_contours(float *hm,int hc,gp_matrix &phi,const char *filename) {
            FILE *fp=gp_safe_fopen(filename,"w");
            output_contours(hm,hc,phi,fp);
            fclose(fp);
        }
        void output_contours(float *hm,int hc,gp_matrix &phi,FILE *fp=stdout);
        void output_bitmap(const char* filename,float zlo,float zhi,int p_ind=0);
        /** Outputs a color bitmap of the field information in PNG format, with
         * one pixel corresponding to each field value. The range of the color
         * paletter is automatically calculated based on the minimum and
         * maximum field values.
         * \param[in] filename the name of the file to write to. */
        inline void output_bitmap(const char* filename,int p_ind=0) {
            float zlo,zhi;
            field_range(zlo,zhi);
            output_bitmap(filename,zlo,zhi,p_ind);
        }
        void output_bitmap_upsample(const char* filename,int ups,bool stagger,float zlo,float zhi,int p_ind=0);
        inline void output_bitmap_upsample(const char* filename,int ups,bool stagger,int p_ind=0) {
            float zlo,zhi;
            field_range(zlo,zhi);
            output_bitmap_upsample(filename,ups,stagger,zlo,zhi,p_ind);
        }
        void field_range(float &zlo,float &zhi);
        void l1_l2_norm(float &l1,float &l2);
        float lp_norm(double p);
    private:
        /** Determines if a grid index is out of bounds.
         * \param[in] (i,j) the grid index to consider.
         * \return Whether the index is out of bounds or not. */
        inline bool oob(int i,int j) {return i<0||j<0||i>=m-1||j>=n-1;}
        int bisection_search(float *q,float t,int lo,int hi);
        bool move(int &i,int &j,bool &vert,bool &back,bool rem);
        void trace_path(int i,int j,bool vert,FILE *fp,float h);
        void trace_path(int i,int j,bool vert,FILE *fp,float h,gp_matrix &phi);
        void ds_row_set(float pre,float *f2,float *f,int m2);
        void ds_row_add(float pre,float *f2,float *f,int m2);
        void ds_halve_range(float *&z,int d,int d2);
        float xsep(int i,int j,float h);
        float ysep(int i,int j,float h);
        inline void general_read(float *ff,int sz,FILE *fp) {
            if(fread(ff,sizeof(float),sz,fp)!=(unsigned int) sz) fatal_error("File input error",3);
        }

        void read_internal(FILE *fp);
        void fatal_error(const char *p,int status);
        void block_indicators(float h);
        void allocate();
        /** Calculates the norm of a given vector.
         * \param[in] (x,y) the vector coordinates.
         * \return The computed norm. */
        inline float dis(float x,float y) {
            return sqrt(x*x+y*y);
        }
        /** Calculates the value of a level set function in another gp_matrix
         * class at an edge position within the contour algorithm.
         * \param[in] (i,j) the grid index associated with the edge.
         * \param[in] vert whether the edge is vertical or not.
         * \param[in] h the value of the contour being plotted. */
        inline float pval(int i,int j,bool vert,float h,gp_matrix &phi) {
            return vert?phi.field_lp(x[i],ysep(i,j,h))
                   :phi.field_lp(xsep(i,j,h),y[j]);
        }
        /** Calculates the position vector on a grid edge associated with a
         * contour.
         * \param[out] (xx,yy) the position vector.
         * \param[in] (i,j) the grid index associated with the edge.
         * \param[in] vert whether the edge is vertical or not.
         * \param[in] h the value of the contour being plotted. */
        inline void pos(double &xx,double &yy,int i,int j,bool vert,float h) {
            if(vert) {xx=x[i];yy=ysep(i,j,h);} else {xx=xsep(i,j,h);yy=y[j];}
        }
        /** Prints the position vector on a grid edge associated with a
         * contour.
         * \param[in] fp the file handle to print to.
         * \param[in] (i,j) the grid index associated with the edge.
         * \param[in] vert whether the edge is vertical or not.
         * \param[in] h the value of the contour being plotted. */
        inline void print_pos(FILE *fp,int i,int j,bool vert,float h) {
            double xx,yy;
            pos(xx,yy,i,j,vert,h);
            fprintf(fp,"%g %g\n",xx,yy);
        }
};

void gp_matrix_bitmap_colorbar(const char* filename,int wid,int hei,int ptype=0,double gamma=1.);

#endif
