#include "gp_matrix.hh"
#include "gp_palette.hh"
#include "gp_pal_defs.hh"

#include <cstring>
#include <limits>
#include <vector>

/** Initializes a gp_matrix class from a file in Gnuplot matrix binary format,
 * setting up the memory to match the file geometry.
 * \param[in] filename the file name to read. */
gp_matrix::gp_matrix(const char *filename) {
    FILE *fp=gp_safe_fopen(filename,"rb");

    // Read in the grid size
    float q[1];
    if(fread(q,sizeof(float),1,fp)!=1) fatal_error("Can't read header",2);
    m=int(*q+0.5);

    // Read x information
    x=new float[m];
    general_read(x,m,fp);

    // Read field information into dynamically extended memory
    std::vector<float*> ft;
    n=0;
    ft.push_back(new float[m+1]);
    while(fread(ft[n],sizeof(float),m+1,fp)==(unsigned int) m+1) {
        ft.push_back(new float[m+1]);n++;
    }
    if(n==0||!feof(fp)) fatal_error("File input error",3);
    fclose(fp);

    // Assemble y and f arrays
    int i,j;
    y=new float[n];
    f=new float[mn=m*n];float *fpp=f;
    for(j=0;j<n;j++) {
        y[j]=*ft[j];
        for(i=1;i<=m;i++) *(fpp++)=ft[j][i];
    }

    // Set spacings
    dx=(x[m-1]-*x)/(m-1);xsp=1/dx;
    dy=(y[n-1]-*y)/(n-1);ysp=1/dy;

    // Clear the temporary memory
    for(i=n;i>=0;i--) delete ft[i];
}

/** Initializes a gp_matrix class, setting up the memory with specific
 * dimensions, but not filling in any data.
 * \param[in] (m_,n_) the grid dimensions. */
gp_matrix::gp_matrix(int m_,int n_) {
    m=m_;n=n_;allocate();
}

/** Initializes a gp_matrix class, setting up the memory with dimensions that
 * match another gp_matrix class.
 * \param[in] gp a reference to the other gp_matrix class. */
gp_matrix::gp_matrix(gp_matrix &gp) {
    m=gp.m;n=gp.n;allocate();
}

/** Initializes a gp_matrix class as a cropped part of another gp_matrix.
 * \param[in] gp the existing gp_matrix to consider.
 * \param[in] (xmin,xmax) the x range to crop to.
 * \param[in] (ymin,ymax) the y range to crop to. */
gp_matrix::gp_matrix(gp_matrix &gp,float xmin,float xmax,float ymin,float ymax) {

    // Find indices
    int ai=bisection_search(gp.x,xmin,0,gp.m-1),
        bi=bisection_search(gp.x,xmax,ai,gp.m-1),
        aj=bisection_search(gp.y,ymin,0,gp.n-1),
        bj=bisection_search(gp.y,ymax,aj,gp.n-1),i,j;

    // Extend to include crossover gridpoints
    ai=ai>0?ai-1:0;bi=bi<gp.m?bi+1:gp.m;
    aj=aj>0?aj-1:0;bj=bj<gp.n?bj+1:gp.n;

    // Compute dimensions and check that they are valid
    m=bi-ai;n=bj-aj;
    if(m<=0||n<=0) fatal_error("Crop dimensions invalid",1);
    allocate();

    // Fill in x and y arrays
    float *p=x,*p2,*p3;
    for(i=ai;i<bi;i++) *(p++)=gp.x[i];
    for(p=y,j=aj;j<bj;j++) *(p++)=gp.y[j];

    // Fill in the data
    for(p=f,j=aj;j<bj;j++) for(p2=gp.f+j*gp.m+ai,p3=p2+m;p2<p3;) *(p++)=*(p2++);
}

/** Sets x and y coordinate arrays to uniformly cover ranges.
 * \param[in] (ax,bx) the x-coordinate range.
 * \param[in] (ay,by) the y-coordinate range. */
void gp_matrix::set_coordinate_ranges(double ax,double bx,double ay,double by) {

    // Set coordinate arrays. Work internally with double precision for higher
    // accuracy, even though the x and y arrays will only be stored in single
    // precision.
    double ddx=(bx-ax)/(m-1),ddy=(by-ay)/(n-1);
    for(int i=0;i<m;i++) x[i]=ax+i*ddx;
    for(int j=0;j<n;j++) y[j]=ay+j*ddy;

    // Set grid spacing constants
    dx=ddx;dy=ddy;xsp=1./ddx;ysp=1./ddy;
}

/** Uses a bisection search to find the position in an ordered array where a
 * given number is located.
 * \param[in] q an array to search through.
 * \param[in] t the value to search for.
 * \param[in] (lo,hi) the lower and upper bounding indexes in the array to
 *                    start the search from.
 * \return The position the given value. */
int gp_matrix::bisection_search(float *q,float t,int lo,int hi) {
    if(q[lo]>=t) return lo;
    if(q[hi]<t) return hi+1;
    int d=hi-lo,e;
    while(d>1) {
        e=lo+(d>>1);
        (q[e]<t?lo:hi)=e;
        d=hi-lo;
    }
    return hi;
}

/** The gp_matrix destructor frees the dynamically allocated memory. */
gp_matrix::~gp_matrix() {
    delete [] f;
    delete [] y;
    delete [] x;
}

/** Allocates the x, y, and field arrays within the class. */
void gp_matrix::allocate() {
    x=new float[m];
    y=new float[n];
    f=new float[mn=m*n];
}

/** Reads in data from a file in Gnuplot matrix binary format. If the
 * horizontal dimension of the file does not match the value of m, then an
 * error is given.
 * \param[in] filename the file name to read. */
void gp_matrix::read(const char *filename) {
    FILE *fp=gp_safe_fopen(filename,"rb");

    // Read in the horizontal grid size and check that it agrees with the
    // current class allocation
    float q[1];
    if(fread(q,sizeof(float),1,fp)!=1) fatal_error("Can't read header",2);
    if(int(*q+0.5)!=m) fatal_error("Size mismatch",6);

    // Read x information
    general_read(x,m,fp);

    // Read in y and field information
    float *yp=y,*fpp=f;
    while(yp<y+n){
        general_read(yp,1,fp);
        general_read(fpp,m,fp);
        yp++;fpp+=m;
    }
    fclose(fp);

    // Set spacings
    dx=(x[m-1]-*x)/(m-1);xsp=1/dx;
    dy=(y[n-1]-*y)/(n-1);ysp=1/dy;
}

/** Transforms the field by taking the power of each value to a given exponent.
 * If the field value is negative, it is set to zero.
 * \param[in] gamma the exponent to use. */
void gp_matrix::power(float gamma) {
    for(float *fp=f,*fe=f+mn;fp<fe;fp++) *fp=*fp<0?-pow(-*fp,gamma):pow(*fp,gamma);
}

/** Crops the field values to a smaller rectangle.
 * \param[in] (ilo,ihi) the horitzontal indices to crop at (lower inclusive,
 *                      upper exclusive).
 * \param[in] (jlo,jhi) the vertical indices to crop at (lower inclusive, upper
 *                      exclusive). */
void gp_matrix::crop(int ilo,int ihi,int jlo,int jhi) {

    // Allocate memory for the cropped field and coordinate ranges
    int m_=ihi-ilo;n=jhi-jlo;mn=m_*n;
    float *x_=new float[m_],*y_=new float[n],*f_=new float[mn];

    // Copy the required part of the coordinate ranges
    memcpy(x_,x+ilo,m_*sizeof(float));
    memcpy(y_,y+jlo,n*sizeof(float));

    // Copy the field values
    for(float *fp=f_,*fs=f+ilo+m*jlo;fp<f_+mn;fp+=m_,fs+=m) memcpy(fp,fs,m_*sizeof(float));

    // Delete memory and update values
    delete [] f;delete [] y;delete [] x;
    m=m_;x=x_;y=y_;f=f_;
}

/** Projects the field values to have zero mean. */
void gp_matrix::zero_project() {
    float sum=0;
    for(float *fp=f,*fe=f+mn;fp<fe;fp++) sum+=*fp;
    sum/=mn;
    for(float *fp=f,*fe=f+mn;fp<fe;fp++) *fp-=sum;
}

/** Outputs the data in the Gnuplot matrix binary format.
 * \param[in] fp the file handle to write to. */
void gp_matrix::output(FILE *fp) {
    float mf[1];*mf=float(m);

    // Output the size and x coordinates
    fwrite(mf,sizeof(float),1,fp);
    fwrite(x,sizeof(float),m,fp);

    // Output the y coordinate and field values
    for(int i=0;i<n;i++) {
        fwrite(y+i,sizeof(float),1,fp);
        fwrite(f+m*i,sizeof(float),m,fp);
    }
}

/** Outputs a double-resolution file where each grid point is converted into a
 * rectangular step.
 * \param[in] fp the file handle to write to. */
void gp_matrix::output_stepped(FILE *fp) {
    int i,j,dm=(m<<1)+1;
    float v,*tmp=new float[dm<n<<1?n<<1:dm],*tp=tmp;

    // Output the size and x coordinates
    *(tp++)=float(m<<1);
    *(tp++)=*x*1.5-x[1]*0.5;
    for(i=0;i<m-1;i++) {
        v=0.5*(x[i]+x[i+1]);
        *(tp++)=v;*(tp++)=v;
    }
    *(tp++)=x[m-1]*1.5-x[m-2]*0.5;
    fwrite(tmp,sizeof(float),dm,fp);

    // Output the y coordinates and field values
    for(j=0;j<n;j++) {
        *tmp=0.5*(j==0?*y*3-y[1]:y[j-1]+y[j]);
        tp=tmp+1;
        for(i=0;i<m;i++) {*(tp++)=f[j*m+i];*(tp++)=f[j*m+i];}
        fwrite(tmp,sizeof(float),dm,fp);
        *tmp=0.5*(j==n-1?3*y[n-1]-y[n-2]:y[j]+y[j+1]);
        fwrite(tmp,sizeof(float),dm,fp);
    }
}

/** Returns a bilinear interpolation of the field at a given location.
 * \param[in] (xi,yi) the position in terms of the integer grid indices at
 *                    which to carry out the interpolation.
 * \return The interpolated value. */
float gp_matrix::field_lp_int(float xi,float yi) {

    // Determine which grid square to consider
    int i=(int) xi,j=(int) yi;
    if(i<0) i=0;else if(i>m-2) i=m-2;
    if(j<0) j=0;else if(j>n-2) j=n-2;
    xi-=i;yi-=j;

    // Carry out the bilinear interpolation in the square
    float *fp=f+(i+m*j);
    return (1-yi)*((1-xi)*(*fp)+xi*fp[1])+yi*((1-xi)*fp[m]+xi*fp[m+1]);
}

/* Calculates an x-transect of the field values.
 * \param[in] ys the y position of the transect.
 * \param[in] xt a pointer in which to store the field values. */
void gp_matrix::x_transect(float ys,float *xt) {

    // Determine which columns of the field values to consider
    ys-=*y;ys*=ysp;
    int j=(int) ys;
    if(j<0) j=0;else if(j>n-2) j=n-2;ys-=j;

    // Carry out linear interpolation of the field values
    float *fp=f+m*j,*xp=xt;
    while(xp<xt+m) {
        *(xp++)=*fp*(1-ys)+fp[m]*ys;
        fp++;
    }
}

/** Calculates a y-transect of the field values.
 * \param[in] xs the x position of the transect.
 * \param[in] yt a pointer in which to store the field values. */
void gp_matrix::y_transect(float xs,float *yt) {

    // Determine which columns of the field values to consider
    xs-=*x;xs*=xsp;
    int i=(int) xs;
    if(i<0) i=0;else if(i>m-2) i=m-2;xs-=i;

    // Carry out linear interpolation of the field values
    float *fp=f+i,*yp=yt;
    while(yp<yt+n) {
        *(yp++)=*fp*(1-xs)+fp[1]*xs;
        fp+=m;
    }
}

/* Calculates the field extrema along a y-transect.
 * \param[in] xs the x position of the transect.
 * \param[out] (fmin,fmax) the minimum and maximum field values. */
void gp_matrix::y_extrema(float xs,float &fmin,float &fmax) {

    // Determine which columns of the field to consider
    xs-=*x;xs*=xsp;
    int i=(int) xs;
    if(i<0) i=0;else if(i>m-2) i=m-2;

    // Determine the minimum and maximum of the field valuesgq
    float *fp=f+i+m,fv;
    fmin=fmax=f[i];
    while(fp<f+mn) {
        fv=*fp*(1-xs)+fp[1]*xs;
        if(fv>fmax) fmax=fv;
        if(fv<fmin) fmin=fv;
        fp+=m;
    }
}

/** Calculates the average field value on a circle.
 * \param[in] (x,y) the center of the circle.
 * \param[in] r the radius of the circle.
 * \return The average value. */
float gp_matrix::circle_average(double x,double y,double r) {
    const int sam=200;
    const double tpi=6.283185307179586476925286766559,h=tpi/sam;
    double th,av=field_lp(x+r,y);
    for(th=h;th<tpi-0.5*h;th+=h) av+=field_lp(x+r*cos(th),y+r*sin(th));
    return av/sam;
}

/** Calculates the average field value on a circle, also storing the minimum and
 * maximum of the sampled values.
 * \param[in] (x,y) the center of the circle.
 * \param[in] r the radius of the circle.
 * \param[out] (mi,mx) the minimum and maximum of the sampled values.
 * \return The average value. */
float gp_matrix::circle_average(double x,double y,double r,double &mi,double &mx) {
    const int sam=200;
    const double tpi=6.283185307179586476925286766559,h=tpi/sam;
    double th,av=field_lp(x+r,y),va;mi=mx=av;
    for(th=h;th<tpi-0.5*h;th+=h) {
        va=field_lp(x+r*cos(th),y+r*sin(th));
        av+=va;
        if(va<mi) mi=va;
        if(va>mx) mx=va;
    }
    return av/sam;
}

/** Calculates the x coordinate of a contour along a horizontal edge.
 * \param[in] (i,j) the grid location to consider.
 * \param[in] h the contour value to consider.
 * \return The x coordinate. */
float gp_matrix::xsep(int i,int j,float h) {
    int ij=i+m*j;
    float df=f[ij]-f[ij+1];
    df=fabs(df)<1e-10?0.5:(f[ij]-h)/df;
    return x[i+1]+(x[i]-x[i+1])*(1-df);
}

/** Calculates the y coordinate of a contour along a vertical edge.
 * \param[in] (i,j) the grid location to consider.
 * \param[in] h the contour value to consider.
 * \return The y coordinate. */
float gp_matrix::ysep(int i,int j,float h) {
    int ij=i+m*j;
    float df=f[ij]-f[ij+m];
    df=fabs(df)<1e-10?0.5:(f[ij]-h)/df;
    return y[j+1]+(y[j]-y[j+1])*(1-df);
}

/** Moves along a contour by one step.
 * \param[in] (i,j) the current grid location.
 * \param[in] vert a boolean value determining whether to start on vertical
 *                 edge.
 * \param[in] back whether to move in the reverse direction.
 * \param[in] rem whether to remove the contour information. */
bool gp_matrix::move(int &i,int &j,bool &vert,bool &back,bool rem) {

    // Loop through the different cases
    if(vert) {
        if(back) {
            i--;
            if(oob(i,j)) {i++;return false;}
            int &v=s[i+m*j];
            switch(v) {
                case 2: vert=false;if(rem) v=0;return true;
                case 3: if(rem) v=0;return true;
                case 7: vert=false;back=false;j++;if(rem) v=0;return true;
                case 6: vert=false;if(rem) v=4;return true;
                case 8: vert=false;back=false;j++;if(rem) v=1;return true;
                default: i++; return false;
            }
        } else {
            if(oob(i,j)) return false;
            int &v=s[i+m*j];
            switch(v) {
                case 1: vert=false;back=true;if(rem) v=0;return true;
                case 3: if(rem) v=0;i++;return true;
                case 4: vert=false;j++;if(rem) v=0;return true;
                case 6: vert=false;j++;if(rem) v=2;return true;
                case 8: vert=false;back=true;if(rem) v=7;return true;
                default: return false;
            }
        }
    } else {
        if(back) {
            j--;
            if(oob(i,j)) {j++;return false;}
            int &v=s[i+m*j];
            switch(v) {
                case 4: vert=true;if(rem) v=0;return true;
                case 5: if(rem) v=0;return true;
                case 6: vert=true;if(rem) v=2;return true;
                case 7: vert=true;back=false;if(rem) v=0;i++;return true;
                case 8: vert=true;back=false;if(rem) v=1;i++;return true;
                default: j++;return false;
            }
        } else {
            if(oob(i,j)) return false;
            int &v=s[i+m*j];
            switch(v) {
                case 1: vert=true;back=true;if(rem) v=0;return true;
                case 2: vert=true;if(rem) v=0;i++;return true;
                case 5: if(rem) v=0;j++;return true;
                case 6: vert=true;if(rem) v=4;i++;return true;
                case 8: vert=true;back=true;if(rem) v=7;return true;
                default: return false;
            }
        }
    }
}

/** Traces out a single contour as a contiguous path.
 * \param[in] (i,j) the starting location to consider.
 * \param[in] vert a boolean value determining whether to start on vertical
 *                 edge.
 * \param[in] fp a file handle to write to.
 * \param[in] h the contour value to consider. */
void gp_matrix::trace_path(int i,int j,bool vert,FILE *fp,float h) {
    int ii=i,jj=j;bool back=true,vvert=vert;

    // Move to one end of the contour (or stop if an entire loop is found)
    while(move(i,j,vert,back,false)) {
        if(i==ii&&j==jj&&vert==vvert) break;}
    back=!back;

    // Trace over the contour and output the positions to a file
    do print_pos(fp,i,j,vert,h);while(move(i,j,vert,back,true));
    fputs("\n\n",fp);
}

/** Traces out a single contour as a contiguous path.
 * \param[in] (i,j) the starting location to consider.
 * \param[in] vert a boolean value determining whether to start on vertical
 *                 edge.
 * \param[in] fp a file handle to write to.
 * \param[in] h the contour value to consider.
 * \param[in] phi a level set function in a gp_matrix class. */
void gp_matrix::trace_path(int i,int j,bool vert,FILE *fp,float h,gp_matrix &phi) {
    int ii=i,jj=j;bool back=true,vvert=vert;

    // Move to one end of the contour (or stop if an entire loop is found)
    while(move(i,j,vert,back,false)) {

        // Search for case when a loop is found
        if(i==ii&&j==jj&&vert==vvert) {

            // If possible, find a point on the loop that is outside the level
            // set, to avoid putting the break in the circular contour at a
            // visible location
            while(move(i,j,vert,back,false)) {
                if(pval(i,j,vert,h,phi)>0) break;
                if(i==ii&&j==jj&&vert==vvert) break;
            }
            break;
        }
    }
    back=!back;

    // Print the first point of the path if it is inside the level set
    ii=i;jj=j;vvert=vert;
    double pp=pval(i,j,vert,h,phi),p=1;
    if(pp<0) print_pos(fp,i,j,vert,h);

    // Trace out the contour
    while(move(i,j,vert,back,true)) {
        p=pval(i,j,vert,h,phi);

        // If the level set switches sign over this segment, then add a point
        // to the file at the switching point
        if((pp<0&&p>=0)||(p<0&&pp>=0)) {
            double x2,y2,x3,y3,z=pp/(pp-p),fz=1-z;
            pos(x2,y2,ii,jj,vvert,h);
            pos(x3,y3,i,j,vert,h);
            fprintf(fp,"%g %g\n",x2*fz+x3*z,y2*fz+y3*z);
            if(pp<0) fputs("\n\n",fp);
        }

        // If the point is inside the level set, then print it
        if(p<0) print_pos(fp,i,j,vert,h);
        ii=i;jj=j;vvert=vert;pp=p;
    }

    // If a contour was being printed, then add a newline to finish it
    if(p<0) fputs("\n\n",fp);
}

/** Outputs the contours of the field, tracing out each contour as a contiguous
 * line for plotting in Gnuplot.
 * \param[in] hm a pointer to an array of contour values to consider.
 * \param[in] hc the number of elements in the array.
 * \param[in] fp the file handle to write to. */
void gp_matrix::output_contours(float *hm,int hc,FILE *fp) {
    s=new int[mn];
    int i,j,l;

    // If multiple fields are requested, then compute the field limits to skip
    // unused cases
    float zlo,zhi;
    if(hc>2) field_range(zlo,zhi);
    else {zhi=std::numeric_limits<float>::max();zlo=-zhi;}

    // Loop over the required contours
    for(float *hp=hm;hp<hm+hc;hp++) {
        float h=*hp;

        // Skip if out of range
        if(h<zlo||h>zhi) continue;

        // Set up block indicators
        block_indicators(h);

        // Output contours
        for(j=0;j<n-1;j++) for(i=0;i<m-1;i++) {
            l=s[i+m*j];
            if(l==1||l==3||l==4||l==6||l==8) trace_path(i,j,true,fp,h);
            if(l==2||l==5) trace_path(i,j,false,fp,h);
            if(l==7) trace_path(i+1,j,true,fp,h);
        }
    }

    // Close file and free dynamically allocated memory
    delete [] s;
}

/** Set up the block indicators array for tracing a contour. For each grid
 * block, store an integer corresponding to what contour segments pass through
 * the block. */
void gp_matrix::block_indicators(float h) {
    int i,j,ij,l;
    for(j=0;j<n-1;j++) for(i=0;i<m-1;i++) {
        ij=i+m*j;
        l=f[ij]>h?1:0;
        if(f[ij+1]>h) l|=2;
        if(f[ij+m]>h) l|=4;
        if(f[ij+m+1]>h) l^=7;
        if(l==6) {
            float x2=xsep(i,j,h),x3=xsep(i,j+1,h);
            float y2=ysep(i,j,h),y3=ysep(i+1,j,h);
            if(dis(x[i]-x2,y[j]-y2)+dis(x[i+1]-x3,y[j+1]-y3)>
               dis(x[i]-x3,y[j]-y3)+dis(x[i+1]-x2,y[j+1]-y2)) l=8;
        }
        s[ij]=l;
    }
}

/** Outputs the contours of the field, tracing out each contour as a contiguous
 * line for plotting in Gnuplot, and removing parts of the contours in regions
 * where a given level set function is positive.
 * \param[in] hm a pointer to an array of contour values to consider.
 * \param[in] hc the number of elements in the array.
 * \param[in] phi a level set function in a gp_matrix class.
 * \param[in] filename the name of the file to write to. */
void gp_matrix::output_contours(float *hm,int hc,gp_matrix &phi,FILE *fp) {
    s=new int[mn];
    int i,j,l;

    // If multiple fields are requested, then compute the field limits to skip
    // unused cases
    float zlo,zhi;
    if(hc>2) field_range(zlo,zhi);
    else {zhi=std::numeric_limits<float>::max();zlo=-zhi;}

    // Loop over the required contours
    for(float *hp=hm;hp<hm+hc;hp++) {
        float h=*hp;

        // Skip if out of range
        if(h<zlo||h>zhi) continue;

        // Set up block indicators
        block_indicators(h);

        // Output contours
        for(j=0;j<n-1;j++) for(i=0;i<m-1;i++) {
            l=s[i+m*j];
            if(l==1||l==3||l==4||l==6||l==8) trace_path(i,j,true,fp,h,phi);
            if(l==2||l==5) trace_path(i,j,false,fp,h,phi);
            if(l==7) trace_path(i+1,j,true,fp,h,phi);
        }
    }

    // Close file and free dynamically allocated memory
    delete [] s;
}

/** Downsamples the field by a factor of two, dealing with cases when the grid
 * has both an odd and even numbe of gridpoints in each dimension. */
void gp_matrix::half_downsample() {
    int j2,m2=(m+1)>>1,n2=(n+1)>>1,mn2=m2*n2;

    // Downsample the two coordinate ranges
    ds_halve_range(x,m,m2);
    ds_halve_range(y,n,n2);

    // Downsample the grid, first taking care of when there an odd number of
    // rows
    float *f2=new float[mn2];
    if(n&1) {
        ds_row_set(1,f2,f,m2);
        for(j2=1;j2<n2-1;j2++) {
            float *fp2=f2+j2*m2,*fp=f+2*j2*m;
            ds_row_set(0.25,fp2,fp-m,m2);
            ds_row_add(0.5,fp2,fp,m2);
            ds_row_add(0.25,fp2,fp+m,m2);
        }
        ds_row_set(1,f2+mn2-m2,f+mn-m,m2);

    // Now deal with the case where there are an even number of rows
    } else for(j2=0;j2<n2;j2++) {
        float *fp2=f2+j2*m2,*fp=f+2*j2*m;
        ds_row_set(0.5,fp2,fp,m2);
        ds_row_add(0.5,fp2,fp+m,m2);
    }

    // Replace the field with the downsampled values and update the grid
    // dimensions
    delete [] f;f=f2;
    m=m2;n=n2;mn=mn2;
    dx*=2;dy*=2;xsp*=0.5;ysp*=0.5;
}

/** Sets a row in the downsampled grid.
 * \param[in] pre a prefactor to apply to the new field values.
 * \param[in] f2 a pointer to the new row.
 * \param[in] f a pointer to the old row.
 * \param[in] m2 the length of the new row. */
void gp_matrix::ds_row_set(float pre,float *f2,float *f,int m2) {
    float *fp=f;
    if(m&1) {
        *f2=pre*(*f);fp+=2;
        for(int i2=1;i2<m2-1;i2++,fp+=2)
            f2[i2]=pre*(0.25*(fp[-1]+fp[1])+0.5*(*fp));
        f2[m2-1]=pre*f[m-1];
    } else for(int i2=0;i2<m2;i2++,fp+=2) f2[i2]=0.5*pre*(*fp+fp[1]);
}

/** Adds to a row in the downsampled grid.
 * \param[in] pre a prefactor to apply to the new field values.
 * \param[in] f2 a pointer to the new row.
 * \param[in] f a pointer to the old row.
 * \param[in] m2 the length of the new row. */
void gp_matrix::ds_row_add(float pre,float *f2,float *f,int m2) {
    float *fp=f;
    if(m&1) {
        *f2+=pre*(*f);fp+=2;
        for(int i2=1;i2<m2-1;i2++,fp+=2)
            f2[i2]+=pre*(0.25*(fp[-1]+fp[1])+0.5*(*fp));
        f2[m2-1]+=pre*f[m-1];
    } else for(int i2=0;i2<m2;i2++,fp+=2) f2[i2]+=0.5*pre*(*fp+fp[1]);
}

/** Cuts a coordinate range in half due to downsampling the field values.
 * \param[in] z a reference to the coordinate range to halve.
 * \param[in] d the current size of the coordinate range.
 * \param[in] d2 the new size of the coordinate range. */
void gp_matrix::ds_halve_range(float *&z,int d,int d2) {
    float *z2=new float[d2];
    if(d&1) {
        for(int i=0;i<d2;i++) z2[i]=z[2*i];
    } else {
        for(int i=0;i<d2;i++) z2[i]=0.5*(z[2*i]+z[2*i+1]);
    }
    delete [] z;z=z2;
}

/** Prints an error message and quits.
 * \param[in] p the message to print.
 * \param[in] status the error code to exit with. */
void gp_matrix::fatal_error(const char *p,int status) {
    fprintf(stderr,"%s\n",p);
    exit(status);
}

/** Outputs a color bitmap of the field information in PNG format, with one
 * pixel corresponding to each field value.
 * \param[in] filename the name of the file to write to.
 * \param[in] (zlo,zhi) the range of field values to use in the color palette.
 * \param[in] p_ind the palette index. */
void gp_matrix::output_bitmap(const char* filename,float zlo,float zhi,int p_ind) {
#ifdef HAS_PNG

    // Create PNG bitmap array
    png_bytep *rowp=new png_bytep[n];

    // Assemble bitmap information
    int i,j;
    float izra=1.0/(zhi-zlo),*fp=f;
    for(j=n-1;j>=0;j--) {
        png_byte* rp=(rowp[j]=new png_byte[3*m]);
        for(i=0;i<m;i++,fp++) gp_pals[p_ind].col((*fp-zlo)*izra,rp);
    }

    // Output the file
    gp_matrix_png_write(m,n,rowp,filename);

    // Free dynamically allocated memory
    for(j=0;j<n;j++) delete [] rowp[j];
    delete [] rowp;
#else

    // Give a status message in the case where the code was compiled without
    // the PNG library
    puts("No PNG support");
#endif
}

/** Outputs a color bitmap of the field information in PNG format, with one
 * pixel corresponding to each field value.
 * \param[in] filename the name of the file to write to.
 * \param[in] ups the upsampling factor to use, so that each field value is
 *                represented by an ups by ups block of pixels.
 * \param[in] stagger whether to render a staggered field, where the field
 *                    values on the boundaries are halved in size; this works
 *                    best when ups is even, so that the halved blocks have a
 *                    clear size.
 * \param[in] (zlo,zhi) the range of field values to use in the color palette.
 * \param[in] p_ind the palette index. */
void gp_matrix::output_bitmap_upsample(const char* filename,int ups,bool stagger,float zlo,float zhi,int p_ind) {
#ifdef HAS_PNG
    printf("%g %g\n",zlo,zhi);

    // Create PNG bitmap array image
    int sta=stagger?-(ups>>1):0,i,j=sta,jj,lj,uj,
        sm=ups*(stagger?m-1:m),sn=ups*(stagger?n-1:n);
    png_bytep *rowp=new png_bytep[sn];

    // Assemble bitmap information
    float izra=1.0/(zhi-zlo),*fp=f;
    for(j=sta;j<sn;) {
        uj=j<=0?sn:sn-j;j+=ups;
        lj=j>=sn?0:sn-j;

        // Allocate memory for one row, and link any duplicate rows to this one
        rowp[lj]=new png_byte[3*sm];
        for(jj=lj+1;jj<uj;jj++) rowp[jj]=rowp[lj];

        // Assemble the bitmap information for one row
        png_byte *rp=rowp[lj],*re;
        for(i=sta;i<sm;fp++) {
            i+=ups;
            re=rowp[lj]+3*(i>=sm?sm:i);
            gp_pals[p_ind].col((*fp-zlo)*izra,rp);
            while(rp<re) {*rp=rp[-3];rp++;}
        }
    }

    // Output the file
    gp_matrix_png_write(sm,sn,rowp,filename);

    // Free dynamically allocated memory
    for(j=sta;j<sn;j+=ups) delete [] rowp[j<0?0:j];
    delete [] rowp;
#else

    // Give a status message in the case where the code was compiled without
    // the PNG library
    puts("No PNG support");
#endif
}

/** Calculates the extramal values of the field.
 * \param[out] zlo the minimum field value.
 * \param[out] zhi the maximum field value. */
void gp_matrix::field_range(float &zlo,float &zhi) {
    zlo=*f;zhi=*f;
    for(float *fp=f+1;fp<f+mn;fp++) {
        if(*fp<zlo) zlo=*fp;
        if(*fp>zhi) zhi=*fp;
    }
}

/** Calculates the L1 and L2 norms of the field.
 * \param[out] l1 the computed L1 norm.
 * \param[out] l2 the computed L2 norm. */
void gp_matrix::l1_l2_norm(float &l1,float &l2) {
    l1=*f;l2=*f*(*f);
    for(float *fp=f+1;fp<f+mn;fp++) {
        l1+=fabs(*fp);l2+=*fp*(*fp);
    }
    float fac=1./mn;
    l1*=fac;l2=sqrt(l2*fac);
}

/** Calculates the Lp norm of the field.
 * \param[in] p the exponent to consider.
 * \return The computed norm. */
float gp_matrix::lp_norm(double p) {
    double lp=pow(fabs(*f),p);
    for(float *fp=f+1;fp<f+mn;fp++) lp+=pow(fabs(*fp),p);
    return static_cast<float>(pow(lp/mn,1./p));
}

/** Calculates the mean value of the field.
 * \return The mean value. */
float gp_matrix::average() {
    float av=*f;
    for(float *fp=f+1;fp<f+mn;fp++) av+=*fp;
    return av/mn;
}

/** Opens a file, and checks the return value to ensure that the operation was
 * successful.
 * \param[in] filename the file to open.
 * \param[in] mode the cstdio fopen mode to use.
 * \return The file handle. */
FILE* gp_safe_fopen(const char *filename,const char *mode) {
    FILE *fp=fopen(filename,mode);
    if(fp==NULL) {
        fprintf(stderr,"Unable to open file '%s'\n",filename);
        exit(1);
    }
    return fp;
}
