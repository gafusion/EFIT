# Makefile for building EFUND
# Revisions:
# $Log: Makefile,v $
# Revision 1.1.2.2  2008/12/10 23:08:33  radhakri
# *** empty log message ***
#
# Revision 1.1.2.1  2008/11/18 22:22:46  radhakri
# *** empty log message ***
#
# Revision 1.1  2007/06/01 05:19:19  renq
#
# Makefile: added by Qilong Ren
#   Makefile for EFUND.
#
# Revision 1.2  2007/05/05 05:28:13  renq
# using make pattern rule
#
# Revision 1.1  2007/04/24 06:03:27  renq
# Initial revision
#
#srini added few lines to include machine dependency
MACHINE := $(shell uname)

#default HP-Unix
FFLAGS = -K +U77 +O1 -w -v -c
LDFLAGS= -K +O1 -w -v -o
FC=f90
FCOM=

ifeq ($(MACHINE),Linux)
#	FC=/c/source/PGI/pgi/linux86/8.0-2/bin/pgf90
#	FC=/c/source/PGI/pgi/linux86/9.0-3/bin/pgf90
#	FC=/c/source/PGI/pgi/linux86/10.8/bin/pgf90
#	FC=/fusion/usc/opt/pgi/linux86-64/15.3/bin/pgf90
	FC=/fusion/usc/opt/pgi/linux86/13.2/bin/pgf90
	FC=/fusion/usc/opt/pgi/linux86-64/15.3/bin/pgf90
	FFLAGS = -Mbyteswapio -fast -c
	LDFLAGS = -Mlfs -o
     ifeq ($(FCOM),lf95)
	FC=/usr/local/lf9562/bin/lf95
	FFLAGS = -c  
	LDFLAGS = -o

     endif
endif
MODULES = exparm.f90 bfgrid.f90 cacoil.f90 cecoil.f90 filech.f90 \
          coilsp.f90 consta.f90 cvesel.f90 fcoil.f90 vas_sub1.f90 \
          fshift.f90 input.f90 mprobe.f90 nio.f90 \
          pmodel.f90 rogowl.f90 siloop.f90 efund_read.f03 efundud.f90
OBJECTS=$(MODULES:.f*=.o)
efundud:$(OBJECTS) 
	$(FC) $(LDFLAGS) $@ $^
	@echo "EFUND BUILT"
%.o:%.f90
	$(FC) $(FFLAGS) $<
clean:
	rm *.o *.mod efundud
	@echo "CLEAN DONE"
