OBJS = chainutil.o quaternion.o inputstruct.o key.o readkey.o genutil.o mt19937.o helix.o energy.o cylinder.o tricubspline.o dbrent.o nrtype.o nrutils.o driver.o optimize.o database.o basinhop.o inputparam.o

FC = gfortran
FLAGS = -O3 
DEBUGFLAGS =
LIBS = -lblas -llapack

.SUFFIXES: .o .f .f90 .mod

.f90.o:
	$(FC) $(FLAGS) $(DEBUGFLAGS) -c $<
.f.o:
	$(FC) $(FLAGS) $(DEBUGFLAGS) -c $<

default: $(OBJS) main.o
	$(FC) $(FLAGS) $(DEBUGFLAGS) $(OBJS) main.o   -o ../fibermodel.exe $(LIBS)

testing: $(OBJS) testmain.o
	$(FC) $(FLAGS) $(DEBUGFLAGS) $(OBJS) testmain.o   -o ../testrun.exe $(LIBS)

clean:
	rm *.o *.mod

chainutil.o: quaternion.o genutil.o key.o 
readkey.o:  genutil.o key.o inputparam.o
inputstruct.o: mt19937.o inputparam.o
helix.o: inputstruct.o tricubspline.o cylinder.o energy.o
genutil.o: mt19937.o
energy.o: chainutil.o quaternion.o cylinder.o tricubspline.o
cylinder.o: tricubspline.o dbrent.o
dbrent.o: nrtype.o nrutils.o
driver.o: chainutil.o cylinder.o energy.o helix.o optimize.o database.o basinhop.o
testmain.o: driver.o
main.o: driver.o
database.o: chainutil.o
basinhop.o: database.o quaternion.o key.o
