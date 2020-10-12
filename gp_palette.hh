#ifndef GP_PALETTE_HH
#define GP_PALETTE_HH

#include <cstdio>
#include <cstring>

#ifdef HAS_PNG
#include <png.h>

void gp_matrix_png_abort(const char* err_msg);
void gp_matrix_png_write(int m_,int n_,png_bytep *rowp,const char *filename);

#endif

class gp_palette {
    public:
        /** Whether the colors are equally spaced. */
        const bool equal;
        /** The number of color reference points. */
        int n;
        /** A pointer to the color information. */
        const float* c;
        /** The short name of the palette. */
        const char* sname;
        /** The long name of the palette. */
        const char* lname;
        gp_palette(bool equal_,int n_,const float *c_,const char* sname_,const char* lname_)
            : equal(equal_), n(n_), c(c_), sname(sname_), lname(lname_) {}
        void output_gnuplot(FILE *fp,double gamma=1);
        void output_text(FILE *fp,double gamma=1);
        void colorbar_horizontal(const char* filename,int wid,int hei,double gamma);
        void colorbar_vertical(const char* filename,int wid,int hei,double gamma);
        inline bool match_short(const char* str) {
            return strcmp(str,sname)==0;
        }
        inline bool match_long(const char* str) {
            return strcmp(str,lname)==0;
        }
        inline float scale_and_clip(float &fv) {
            fv=(1/255.)*(fv*256-0.5);
            return fv<0?0:(fv>1?1:fv);
        }
        /** Computes an (R,G,B) color from a scalar input.
         * \param[in] fv the scalar.
         * \param[out] (re,gr,bl) the components of the color. */
        inline void col(float fv,float &re,float &gr,float &bl) {
            col_internal(scale_and_clip(fv),re,gr,bl);
        }
#ifdef HAS_PNG
        /** Computes a PNG color from a scalar input.
         * \param[in] fv the scalar.
         * \param[in,out] rp a pointer to where to store the PNG color, incremented to
         *                   the next entry upon completion. */
        inline void col(float fv,png_byte*& rp) {
            float re,gr,bl;
            col(scale_and_clip(fv),re,gr,bl);
            *(rp++)=static_cast<png_byte>(re+0.5);
            *(rp++)=static_cast<png_byte>(gr+0.5);
            *(rp++)=static_cast<png_byte>(bl+0.5);
        }
#endif
    private:
        void col_internal(float fv,float &re,float &gr,float &bl);
        inline unsigned int clip(float fv) {
            int v=static_cast<int>(fv+0.5);
            return v<0?0:(v>255?255:v);
        }
        inline unsigned int hex_col(const float *cp) {
            return clip(cp[2])|(clip(cp[1])<<8)|(clip(*cp)<<16);
        }
};

#endif
