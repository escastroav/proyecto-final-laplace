CXX = g++
OTHFLAGS = -O0 -Wall 
CXXFLAGS = -I $HOME/local/include -Wl,-rpath=$HOME/local/lib -std=c++11 
LDFLAGS =  -L $HOME/local/lib -lpapi 
SOURCES = laplace-cz.cpp

OBJ = $(SOURCES:.cpp=.o) 
DEPS = 

all : main.x $(SOURCES) $(DEPS) # target by default

main.x : $(OBJ)
	@echo "Creating THE main MULTIPLICATION executable..."
	$(CXX) $(OTHFLAGS) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

.cpp.o:
	$(CXX) -c -o $@ $^ $(CXXFLAGS)

.PHONY: clean
clean: 
	rm -f *.o *~ *.x
