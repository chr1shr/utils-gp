#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>

#include "gp_pal_defs.hh"
#include "gp_palette.hh"

void syntax_message() {
    fputs("Usage: ./gpp_info list\n"
          "       ./gpp_info gnu <palette>\n"
          "       ./gpp_info txt <palette>\n",stderr);
    exit(1);
}

int main(int argc,char **argv) {

	// Check for at least one command-line argument
	if(argc<2) syntax_message();

    // Check for the command to list the palettes
    if(strcmp(argv[1],"list")==0) {
        if(argc!=2) syntax_message();
        gp_list_palettes(stdout);
        return 1;
    }

    // Check the either Gnuplot or text output is requested
    if(argc!=3) syntax_message();
    bool gnu;
    if(strcmp(argv[1],"gnu")==0) gnu=true;
    else if(strcmp(argv[1],"txt")==0) gnu=false;
    else syntax_message();

	// Find the palette type
    int p_ind=gp_pal_index(argv[2]);
    if(p_ind==-1) {
        fputs("Unknown palette type; available palettes are:\n\n",stderr);
        gp_list_palettes(stderr);
        return 1;
    }

    // Output the palette in the requested format
    gnu?gp_pals[p_ind].output_gnuplot(stdout)
       :gp_pals[p_ind].output_text(stdout);
}
