#!/usr/bin/perl
open A,"palettes.txt" or die "Can't open palettes file\n";

$p=0;
@c=();
while(<A>) {
    $l++;

    # Skip any blank lines
    next if /^[ \t]*\n$/;

    # Search for a palette header line
    if(/# (.*) \[(...)\]/) {

        # Store palette header information
        $lf[$p]=$1;$sf[$p]=$2;
        $ke[$p]=$#c+1;$nu[$p]=1;

        # Process the first color. If there are three entries, this palette has
        # equal spacing between colors; if there are four entries this palette
        # has unequal spacing.
        @a=split ' ',<A>;$l++;
        die "Can't parse color on line $l\n" unless $#a==2||$#a==3;
        $n[$p]=$#a;
        push @c,@a;

        # Search for additional colors
        while(<A>) {
            $l++;
            
            # Check for completion of this palette
            last if /^[ \t]*\n$/;

            # Parse the color, and make sure that the number of entries matches
            # the first color
            @a=split;$nu[$p]++;
            die "Can't parse color on line $l\n" unless $#a==2||$#a==3;
            die "Inconsistent number of entries on line $\n" unless $#a=$n[$p];
            push @c,@a;
            die "$@ at line $l\n" if $@;
        }
        $p++;
    } else {
        die "Can't parse header on line $l\n";
    }
}
close A;

# Output the C++ header file, declaring the two palette data structures
$q=$#c+1;
open A,">gp_pal_defs.hh" or die "Can't open palette info header file\n";
print A <<EOF;
#ifndef GP_PAL_DEFS_HH
#define GP_PAL_DEFS_HH

#include <cstdio>

#include "gp_palette.hh"

extern const float gp_pal_cols[$q];
extern gp_palette gp_pals[$p];

int gp_pal_index(const char *str);
void gp_list_palettes(FILE *fp=stdout);

#endif
EOF
close A;

# Output the C++ implementation file
open A,">gp_pal_defs.cc" or die "Can't open palette info implementation file\n";
print A <<EOF;
#include "gp_pal_defs.hh"

const float gp_pal_cols[$q]={
EOF

# Output the palette colors in a compressed format
$k=0;
while($k<$q-12) {
    print A "    ";
    print A "$c[$_]," foreach $k..($k+10);
    $k+=11;
    print A "$c[$k],\n";
    $k++;
}
print A "    ";
print A "$c[$_]," foreach $k..($q-2);
$k=$q-1;

# Output the palette descriptions
print A "$c[$k]\n};\n\ngp_palette gp_pals[$p]={\n";
printf A '    gp_palette(%s,%d,gp_pal_cols+%d,"%s","%s")%s',
    $n[$_]==2?"true":"false",$nu[$_],$ke[$_],
    $lf[$_],$sf[$_],$_==$p-1?"\n":",\n" foreach (0..$p-1);

# Output the function for matching palette names 
print A <<EOF;
};

/** Finds the index of the palette matching a given string.
 * \\param[in] the string to match.
 * \\return The palette index; if no palette matches then -1 is returned. */
int gp_pal_index(const char *str) {
    for(int i=0;i<$p;i++)
        if(gp_pals[i].match_short(str)) return i;
    for(int i=0;i<$p;i++)
        if(gp_pals[i].match_long(str)) return i;
    return -1;
}

/** Lists the available palettes.
 * \\param[in] fp the file handle to write to. */
void gp_list_palettes(FILE *fp) {
    for(int i=0;i<$p;i++)
        fprintf(fp,"%d: %s [%s]\\n",i,gp_pals[i].sname,gp_pals[i].lname);
}
EOF

close A;
