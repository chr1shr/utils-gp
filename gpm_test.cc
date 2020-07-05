#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>

#include "gp_matrix.hh"

int main() {
    const int m=11,n=11;
    gp_matrix gp(m,n);

    for(int j=0;j<n;j++) {
        gp.y[j]=0.1*j;
        for(int i=0;i<m;i++) {
            gp.x[i]=0.1*i;
            gp.f[i+j*m]=0.1*i;
        }
    }

    gp.output("f.0");
    gp.half_downsample();
    gp.output("fh.0");
}
