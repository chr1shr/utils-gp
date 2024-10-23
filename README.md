# utils-gp: A collection of utilities for analyzing Gnuplot matrix files
[Gnuplot](http://www.gnuplot.info) is a freeware graphing utility that is
widely used and available on many computing platforms. Gnuplot can work with
many types of input data in both text and binary formats.

For plotting two-dimensional fields on a rectangular grid, Gnuplot uses a
"binary matrix" format, which is a compact and efficient way to store both the
field values and the coordinate ranges. This Git repository contains a variety
of utilities for working with data files in this format.

The utilities are structured around a C++ class **gp_matrix** for reading in
these data files and performing a variety of analyses. This class can be used
and linked from other programs. A key feature of the class is that it can deal
with [level set functions](https://en.wikipedia.org/wiki/Level-set_method) that
define one-dimensional interfaces, and it contains several routines that are
tailored toward typical problems that emerge when using level set methods.

Several executables are provided, which can

- make pixel-accurate bitmaps using a variety of custom palettes;
- clean fields based on clipping with a level set function;
- trace contours of the field, including contours that are clipped by a level
  set function.

# The binary matrix format
The binary matrix format is described by typing "help binary matrix
nonuniform" within Gnuplot. Consider a two-dimensional m &times; n grid with
field values f<sub>i,j</sub> for i=0,...,m-1 and j=0,...,n-1. Let
x<sub>i</sub> be the x coordinates of the grid, and y<sub>j</sub> be the
y coordinates of the grid. The format consists of single-precision (32 bit)
floating point numbers, arranged in the following way:

|m|x<sub>0</sub>|x<sub>1</sub>|...|x<sub>m-1</sub>|
|-|-------------|-------------|---|---------------|
|y<sub>0</sub>|f<sub>0,0</sub>|f<sub>1,0</sub>|...|f<sub>m-1,0</sub>|
|y<sub>1</sub>|f<sub>0,1</sub>|f<sub>1,1</sub>|...|f<sub>m-1,1</sub>|
|...|...|...| |...|
|y<sub>n-1</sub>|f<sub>0,n-1</sub>|f<sub>1,n-1</sub>|...|f<sub>m-1,n-1</sub>|

The very first value contains m, the width of the grid, converted to a
floating point number. After reading this value, Gnuplot knows to expect that
the following m entries are x coordinates. Gnuplot can then process batches of
(m+1) entries, as a y coordinate followed by a row of field values. The value
of n is determined once the end of the data file is reached.

While most scientific simulations internally work with double-precision (64
bit) floating numbers, single-precision is usually adequate for plotting and
analysis, and cuts down the data storage requirements by a factor of two. Since
single precision floating point numbers are accurate to eight significant
figures, the rounding errors have an indistinguishable effect on plotted lines.

# Compiling the code
The code is written in C++ and has been tested on Linux, MacOS, and Windows
(via [Cygwin](https://www.cygwin.com)). It uses [Perl](http://www.perl.org)
during the compilation procedure. The following documentation assumes you are
familiar with the Linux/Mac/Cygwin
[command-line interface](https://en.wikipedia.org/wiki/Command-line_interface).

The code requires [libpng](http://www.libpng.org/pub/png/) for making PNG
bitmaps. However, this dependency is not necessary to compile the code. libpng
is available as a standard downloadable package on most systems.

To compile the code it is necessary to create a common configuration file
called **config.mk** in the parent directory, which can be used by all three
repositories. Several templates are provided in the **config** directory. To
use, copy one of the templates into the parent directory. From the utils-gp
directory, on a Linux computer, type
```Shell
cp config/config.mk.linux ../config.mk
```
On a Mac using GCC 14 installed via [MacPorts](http://www.macports.org), type
```Shell
cp config/config.mk.mac_mp ../config.mk
```
On a Mac using GCC installed via [Homebrew](http://brew.sh), type
```Shell
cp config/config.mk.mac_hb ../config.mk
```
On a Windows computer with Cygwin installed, type
```Shell
cp config/config.mk.win_cw ../config.mk
```
In all of these configuration files, PNG support is enabled with the
"-DHAS\_PNG" compiler flag. If libpng is unavailable then this flag can be
removed. After this, the code can be compiled by typing
```Shell
make
```
This will build the a static library called **libgpmtx.a** that can be linked
to by other programs. It will also build the executables **bitmap_field**,
**colorbar**, **gen_test**, **gpm_extrema**, **gpm_metrics**, **gpm_process**,
**gpp_info**, **make_contour**, and **output_text**. Running each executable
gives a short message about syntax.

## Examples
### Generating some test files
A small program called **gen_test** is provided that creates some samples
files in the Gnuplot matrix binary format. This can be run with command
```Shell
./gen_test 100
```
This will create two fields on a 100 &times; 100 grid with coordinate ranges 0
&le; x &le; 1 and 0 &le; y &le; 1:

- waves.fe which contains a field pattern of waves
- sellipse.fe, whose zero contour is a superellipse, (x-0.5)<sup>4</sup> +
  (y-0.5)<sup>4</sup> = 0.4<sup>4</sup>.

The waves field pattern can be examined the Gnuplot command prompt by typing
```Gnuplot
# 3D wireframe plot
set view 30,340
set hidden3d
set xlabel 'x'
set ylabel 'y'
set zlabel 'z'
splot 'waves.fe' matrix binary with lines

# Heat map plot
set pm3d map
plot [0:1] [0:1] 'waves.fe' matrix binary with image
```
These two images are shown below.

![3D plot (left) and heat map plot (right) of the waves test file](http://math.lbl.gov/~chr/utils-gp/gt1.png)

### Palettes
The tools contains a number of built-in palettes for plotting two-dimensional
fields. A utility **gpp_info** is provided for obtaining data about the
palettes. The following command will list the available palettes:
```Shell
./gpp_info list
```
Each palette can be referenced by its three-letter code or full name. A
horizontal color bar for the "heat" palette can be created with the command
```Shell
./color_bar h heat 600 100 heat.png
```

### Creating bitmaps
A test field with higher resolution can be created with the command
```Shell
./gen_test 400
```
The following command will create a bitmap of the waves test field:
```Shell
./bitmap_field -p hea waves.fe waves.png
```
Each pixel in this image corresponds to exactly a single grid point in the
input field. By default the utility stretches the color range to exactly match
the range of the input data. However, two additional arguments can be provided
to the utility to provide the minimum and maximum of the range.

In addition the waves test field can be clipped by the superellipse field
using the command
```Shell
./gpm_process -c sellipse.fe waves.fe waves_clip.fe
```
After this, the clipped field can be plotted in an alternative color
scheme using
```Shell
./bitmap_field -p psy waves_clip.fe waves_clip.png
```
The two resulting images are shown below.

![Two bitmap images made from the waves test file](http://math.lbl.gov/~chr/utils-gp/gt2.png)

### Making a contour
Contours of the waves test field can be made with the command
```Shell
./make_contour waves.fe waves.ctr -1 -0.5 0 0.5 1
```
Here, five contours will be traced out at -1, -0.5, 0, 0.5, and 1. If "r" is
detected as a command-line argument, the utility accepts an alternative syntax
providing (a) the starting contour value, (b) the contour increment, and (c)
the number of contours. Hence the following command is equivalent:
```Shell
./make_contour waves.fe waves.ctr r -1 0.5 5
```
The contours can be plotted in Gnuplot using the commands
```Gnuplot
set xlabel 'x'
set ylabel 'y'
plot [0:1] [0:1] 'waves.ctr' with lines
```
The superellipse zero contour can be calculated using
```Shell
./make_contour sellipse.fe sellipse.ctr
```
Note that if no contour values are given, then the utility just computes the
zero contour. The waves contour can be trimmed to the superellipse
with the command
```Shell
./make_contour -p sellipse.fe waves.fe waves_trim.ctr
```
The trimmed contours and the superellipse can be plotted using the command
```Gnuplot
plot [0:1] [0:1] 'waves_trim.ctr' with lines, 'sellipse.ctr' with lines
```
Plots of the untrimmed and trimmed contours are shown below.

![Original contours (left) and trimmed contours (right) of the waves test file](http://math.lbl.gov/~chr/utils-gp/gt3.png)


## Contact
For questions about the code, contact [Chris Rycroft](http://seas.harvard.edu/~chr/).

## References that use these utilities
This set of utilities has been developed by Chris Rycroft since 2012 and have
been used in a variety of research projects featuring two-dimensional fields.
The utilities are designed to supplement the plotting capabilities of Gnuplot.
They are designed to work rapidly, and are primarily aimed for use in scripting
where they can be used to batch-process many files (_e.g._ to make frames for a
movie).

The following academic articles have made use of these utilities:

1. Chris H. Rycroft and Frédéric Gibou, *Simulations of a stretching bar using
   a plasticity model from the shear transformation zone theory*, J. Comput.
   Phys. **231**, 2155–2179 (2012).
   [doi:10.1016/j.jcp.2011.10.009](https://doi.org/10.1016/j.jcp.2011.10.009)

2. Ken Kamrin, Chris H. Rycroft, and Jean-Christophe Nave, *Reference map
   technique for finite-strain elasticity and fluid–solid interaction*, J.
   Mech. Phys. Solids. **60**, 1952–1969 (2012).
   [doi:10.1016/j.jmps.2012.06.003](https://doi.org/10.1016/j.jmps.2012.06.003)

3. Boris Valkov, Chris H. Rycroft, and Ken Kamrin, *Eulerian method for
   multiphase interactions of soft solid bodies in fluids*, J. Appl. Mech.
   **82**, 041011 (2015).
   [doi:10.1115/1.4029765](https://doi.org/10.1115/1.4029765)

4. Chris H. Rycroft, Yi Sui, and Eran Bouchbinder, *An Eulerian projection
   method for quasi-static elastoplasticity*, J. Comput. Phys. 300, 136–166
   (2015).
   [doi:10.1016/j.jcp.2015.06.046](https://doi.org/10.1016/j.jcp.2015.06.046)

5. Chris H. Rycroft and Martin Z. Bazant, *Asymmetric collapse by dissolution
   or melting in a uniform flow*, Proc. Roy. Soc. A **472**, 20150531 (2016).
   [doi:10.1098/rspa.2015.0531](https://doi.org/10.1098/rspa.2015.0531)

6. Manish Vasoya, Chris H. Rycroft, and Eran Bouchbinder, *Notch fracture
   toughness of glasses: Rate, age and geometry dependence*, Phys. Rev. Applied
   **6**, 024008 (2016).
   [doi:10.1103/PhysRevApplied.6.024008](https://doi.org/10.1103/PhysRevApplied.6.024008)

7. Christoph A. Weber, C. H. Rycroft, and L. Mahadevan, *Differential Activity
   drives Instabilities in Biphasic Active Matter*, Phys. Rev. Lett. **120**,
   248003 (2018).
   [doi:10.1103/PhysRevLett.120.248003](https://doi.org/10.1103/PhysRevLett.120.248003)

8. Chris H. Rycroft, Chen-Hung Wu, Yue Yu, and Ken Kamrin, *Reference map
   technique for incompressible fluid–structure interaction*, J. Fluid Mech.
   **898**, A9 (2020).
   [doi:10.1017/jfm.2020.353](https://doi.org/10.1017/jfm.2020.353)
