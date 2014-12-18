#------ GNU g++ -----------------------------------------------------

CC         = g++
CPPFLAGS     = -Wall -O3 -I${CODEDIR} -Iinclude

LINKER      = g++
LDFLAGS     = -lboost_program_options -lgsl -lgslcblas

#--------------------------------------------------------------------

MAKEDEPEND  = makedepend
DEPENDFLAGS = -m 

#--------------------------------------------------------------------

SOURCES     = solve_ode.cc
OBJECTS     = $(patsubst %.cc,out/%.o,${SOURCES})

EXEC        = solve_ode.x
#--------------------------------------------------------------------

all: bin/$(EXEC)

bin/$(EXEC): $(OBJECTS) 
	@echo ... linking $<: 
	@echo $(OBJECTS)
	$(LINKER) $(OBJECTS) $(LDFLAGS) $(OPTFLAGS) $(CPPFLAGS) -o $@
	@echo

out/%.o: src/%.cc include/%.hh
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

out/solve_ode.o: include/file_io.hh include/tryps.hh include/pa_macros.hh \
                 include/ode_solver.hh include/ode_io_utils.hh include/convergence.hh
