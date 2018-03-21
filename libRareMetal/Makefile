SUBDIRS=mvtnorm pdf src

PARENT_MAKE := Makefile.tool
include Makefile.inc

pdf: src

clean:$(SUBDIRS)
	rm -f $(PROG_OPT)
	rm -f $(PROG_DEBUG)
	rm -f $(PROG_PROFILE)
