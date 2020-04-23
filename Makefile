CXX=mpicxx
MPIF90=mpiifort
CXXFLAGS+=-I./include
CXXFLAGS+=-std=c++11
LIBPATH+=
vpath %.cpp src
vpath %.h include
vpath %.o obj
INJECTION: main.o constant.o indexref.o readqebands.o banddot.o
	mkdir -p obj bin
	$(CXX) -o injection $(LIBPATH) main.o constant.o indexref.o readqebands.o banddot.o
	mv *o bin/
%.o:%.cpp $(DEPS)
	$(CXX) $(CXXFLAGS) -c -o $@ $^
%.o:%.f90 $(DEPS)
	$(MPIF90) -c -o $@ $^
clean:
	rm -rf bin obj *.o
