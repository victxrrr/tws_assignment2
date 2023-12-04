FC=gfortran
EXT=f90
SRC=src/

ifeq ($(FC),gfortran)
	MOD_FLAG_OUT := -J$(SRC)
else ifeq ($(FC),ifort)
	MOD_FLAG_OUT := -module $(SRC)
else ifeq ($(FC),nagfor)
    MOD_FLAG_OUT := -mdir $(SRC)
else ifeq ($(FC),g95)
	MOD_FLAG_OUT := -fmod=$(SRC)
endif

ifeq ($(FC),gfortran)
	MOD_FLAG_IN := -I$(SRC)
else ifeq ($(FC),ifort)
	MOD_FLAG_IN := -I$(SRC)
else ifeq ($(FC),nagfor)
    MOD_FLAG_IN := -I $(SRC)
else ifeq ($(FC),g95)
	MOD_FLAG_IN := -I$(SRC)
endif


FFLAGS_g95      = -g -pedantic -Wall -fbounds-check -ftrace=full
FFLAGS_gfortran = -g -pedantic -Wall -Wunderflow -fbounds-check -fimplicit-none 
FFLAGS_ifort    = -g -debug full -implicitnone -check -free
FFLAGS_nagfor   = -g -gline -u -info -colour

# Select the right flags for the current compiler
FFLAGS=$(FFLAGS_$(FC))
# Comment out to use optimization flags
# FFLAGS += -O2

.PHONY: plot

all: q3

utils.o: $(SRC)utils.$(EXT)
	$(FC) $(FFLAGS) $(MOD_FLAG_OUT) -c $^ 
						
q2.o: $(SRC)q2.$(EXT)
	$(FC) $(FFLAGS) $(MOD_FLAG_IN) -c $^

q2: utils.o q2.o
	$(FC) $(FFLAGS) $(MOD_FLAG_IN) -o $@.out $^ -lblas -llapack

tweaking.o: $(SRC)tweaking.$(EXT)
	$(FC) $(FFLAGS) $(MOD_FLAG_IN) -c $^

q3: utils.o tweaking.o 
	$(FC) $(FFLAGS) $(MOD_FLAG_IN) -o $@.out $^ -lblas -llapack

stability.o: $(SRC)stability.$(EXT)
	$(FC) $(FFLAGS) -c $^
	
q4: stability.o
	$(FC) $(FFLAGS) -o $@.out $^ -lblas -llapack

# to avoid errors when using solver_gfortran.o
%.o: %.mod

plot:
	cd ./plot && pdflatex plot.tex >/dev/null && cd ..

clean:
	@ rm -f *.o *.out *.mod src/*.mod