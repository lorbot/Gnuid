# $Id: Makefile 3101 2008-10-16 19:34:03Z roystgnr $


# The location of the mesh library
LIBMESH_DIR = ~/src2/libmesh-0.7.3.1/libmesh/

# include the library options determined by configure.  This will
# set the variables INCLUDE and LIBS that we will need to build and
# link with the library.
include $(LIBMESH_DIR)/Make.common


###############################################################################
# File management.  This is where the source, header, and object files are
# defined

#
# source files
srcfiles 	:= $(wildcard *.C)

#
# object files
objects		:= $(patsubst %.C, %.$(obj-suffix), $(srcfiles))
###############################################################################



.PHONY: clean clobber distclean

###############################################################################
# Target:
#
target 	   := ./gnuid-$(METHOD)


all:: $(target)

# Production rules:  how to make the target
$(target): $(objects)
	@echo "Linking "$@"..."
	@$(libmesh_CXX) $(libmesh_CXXFLAGS) $(objects) -o $@ $(libmesh_LIBS) $(libmesh_LDFLAGS)


# Useful rules.
clean:
	@rm -f $(objects) *~

clobber:
	@$(MAKE) clean
	@rm -f $(target)

distclean:
	@$(MAKE) clobber
	@rm -f *.o *.g.o *.pg.o

run: $(target)
	@echo "***************************************************************"
	@echo "* Running Example " $(LIBMESH_RUN) $(target) -d 3 $(LIBMESH_DIR)/reference_elements/3D/one_hex27.xda $(LIBMESH_OPTIONS)
	@echo "***************************************************************"
	@echo " "
	@$(LIBMESH_RUN) $(target) -d 3 $(LIBMESH_DIR)/reference_elements/3D/one_hex27.xda $(LIBMESH_OPTIONS)
	@echo " "
	@echo "***************************************************************"
	@echo "* Done Running Example " $(LIBMESH_RUN) $(target) -d 3 $(LIBMESH_DIR)/reference_elements/3D/one_hex27.xda $(LIBMESH_OPTIONS)
	@echo "***************************************************************"


# include the dependency list
include .depend


#
# Dependencies
#
.depend:
	@$(perl) $(LIBMESH_DIR)/contrib/bin/make_dependencies.pl -I. $(foreach i, $(wildcard $(LIBMESH_DIR)/include/*), -I$(i)) "-S\$$(obj-suffix)" $(srcfiles) > .depend
	@$(perl) -pi -e 's#    $(LIBMESH_DIR)#    \$$\(LIBMESH_DIR\)#' .depend
	@echo "Updated .depend"

###############################################################################
