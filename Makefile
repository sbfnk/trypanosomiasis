#------ GNU g++ -----------------------------------------------------

CC         = g++
CPPFLAGS     = -Wall -Iinclude -g # -O3

LINKER      = g++
LDFLAGS     = -lboost_program_options -lgsl -llapack -g

#--------------------------------------------------------------------

MAKEDEPEND  = makedepend
DEPENDFLAGS = -m 

#--------------------------------------------------------------------

SOURCES     = reservoirs.cc ihs.cc
OBJECTS     = $(patsubst %.cc,out/%.o,${SOURCES})

EXEC        = reservoirs
#--------------------------------------------------------------------

all:  bin/$(EXEC)

bin/$(EXEC): $(OBJECTS) 
	@echo ... linking:
	$(LINKER) $(OBJECTS) $(LDFLAGS) $(OPTFLAGS) $(CPPFLAGS) -o $@
	@echo

out/%.o : src/%.cc include/%.hh
	@echo ... compiling $<:
	$(CC) -c $(CPPFLAGS) $(OPTFLAGS) $< -o $@
	@echo

depend:
	$(MAKEDEPEND) $(DEPENDFLAGS) src/*.cc

htmldoc:
	doxygen

pdfdoc: htmldoc
	$(MAKE) -C doc/latex pdf

clean:
	@echo ... cleaning $(OBJECTS) $(EXEC)
	@rm -f bin/* out/*
	@echo ... done

cleaner:
	@echo '... cleaning also *~'
	@rm -f bin/* out/* *~ .#*
	@echo ... done

#-------------------------------------------------------------------
# DO NOT DELETE
