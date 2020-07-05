include ../config.mk

# Lists of files to be built
objs=gp_matrix.o gp_palette.o gp_pal_defs.o
src=$(patsubst %.o,%.cc,$(objs))
execs=gen_test gpm_process gpm_metrics make_contour bitmap_field \
	  colorbar gpp_info gpm_extrema output_text

all: lib $(execs)

include Makefile.dep

depend:
	$(cxx) -MM $(src) >Makefile.dep

clean:
	rm $(objs) $(execs)

gp_pal_defs%hh gp_pal_defs%cc: palettes.txt pal_parse.pl
	perl pal_parse.pl

lib: libgpmtx.a

libgpmtx.a: $(objs)
	rm -f libgpmtx.a
	ar rs libgpmtx.a $^

output_text: output_text.cc libgpmtx.a
	$(cxx) $(cflags) -o $@ $< -L. -lgpmtx $(png_lflags)

gen_test: gen_test.cc libgpmtx.a
	$(cxx) $(cflags) -o $@ $< -L. -lgpmtx $(png_lflags)

gpm_process: gpm_process.cc libgpmtx.a
	$(cxx) $(cflags) -o $@ $< -L. -lgpmtx $(png_lflags)

gpm_metrics: gpm_metrics.cc libgpmtx.a
	$(cxx) $(cflags) -o $@ $< -L. -lgpmtx $(png_lflags)

gpp_info: gpp_info.cc libgpmtx.a
	$(cxx) $(cflags) -o $@ $< -L. -lgpmtx $(png_lflags)

gpm_extrema: gpm_extrema.cc libgpmtx.a
	$(cxx) $(cflags) -o $@ $< -L. -lgpmtx $(png_lflags)

make_contour: make_contour.cc libgpmtx.a
	$(cxx) $(cflags) -o $@ $< -L. -lgpmtx $(png_lflags)

bitmap_field: bitmap_field.cc libgpmtx.a
	$(cxx) $(cflags) -o $@ $< -L. -lgpmtx $(png_lflags)

colorbar: colorbar.cc libgpmtx.a
	$(cxx) $(cflags) -o $@ $< -L. -lgpmtx $(png_lflags)

%.o: %.cc
	$(cxx) $(cflags) $(png_iflags) -c $<

.PHONY: clean depend lib
