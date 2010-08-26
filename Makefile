#------ GNU g++ -----------------------------------------------------

CC         = g++
CPPFLAGS     = -Wall -g -O3 -Iinclude

LINKER      = g++
LDFLAGS   = 

#--------------------------------------------------------------------

MAKEDEPEND  = makedepend
DEPENDFLAGS = -m 

#--------------------------------------------------------------------

SOURCES     = tryps.cc 
OBJECTS     = $(patsubst %.cc,out/%.o,${SOURCES})

EXEC        = tryps
#--------------------------------------------------------------------

all:  bin/$(EXEC)

bin/$(EXEC): $(OBJECTS) 
	@echo ... linking:
	$(LINKER) $(LDFLAGS) $(OPTFLAGS) $(CPPFLAGS) $(OBJECTS) -o $@
	@echo

out/%.o : src/%.cc include/%.hh
	@echo ... compiling $<:
	$(CC) -c $(CPPFLAGS) $(OPTFLAGS) $< -o $@
	@echo

depend:
	$(MAKEDEPEND) $(DEPENDFLAGS) $(SOURCES)

clean:
	@echo ... cleaning $(OBJECTS) $(EXEC)
	@rm -f bin/* out/*
	@echo ... done

cleaner:
	@echo '... cleaning also *~, *.x and gmon.out
	@rm -f bin/* out/* *~ *.x gmon.out $(EXEC)
	@echo ... done

#-------------------------------------------------------------------
# DO NOT DELETE

tryps.o: tryps.hh
