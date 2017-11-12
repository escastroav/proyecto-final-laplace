CXX = g++
OTHFLAGS = -O0 -Wall -fopenmp
CXXFLAGS = -I $HOME/local/include -Wl,-rpath=$HOME/local/lib -std=c++11 
LDFLAGS =  -L $HOME/local/lib -lpapi 
SOURCES = main.cpp laplace.cpp

OBJ = $(SOURCES:.cpp=.o) 
DEPS = laplace.h

all : main.x $(SOURCES) $(DEPS) # target by default

main.x : $(OBJ)
	@echo "Creating THE main executable..."
	$(CXX) $(OTHFLAGS) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

.cpp.o:
	$(CXX) -c -o $@ $^ $(CXXFLAGS)

.PHONY: clean
clean: 
	rm -f *.o *~ *.x
