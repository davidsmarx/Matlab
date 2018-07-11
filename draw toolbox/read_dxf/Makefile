# This Makefile has been tested on Linux only
# but, it should work perfectly well on other systems
TARGET = read_dxf.mexglx

SOURCES = read_dxf.cpp dxf.cpp dxftoken.cpp

# The MATLAB base directory
MATDIR = /usr/local/matlab14
INCDIR = -I. -I$(MATDIR)/extern/include

# the followign line may need to be changed
# for older MATLAB versions or other compile 
# systems
LIBDIR = -L$(MATDIR)/bin/glnx86
LIBS = -lm -lmx -lmat

# The C++ compiler and flags
CXX = g++
CXXFLAGS = -O3 -fomit-frame-pointer

$(TARGET) : $(patsubst %.cpp,%.o,$(SOURCES))
	$(CXX) $(CXXFLAGS) $(LIBDIR) $(LIBS) -shared $^ -o $@

%.o : %.cpp %.h
	$(CXX) $(CXXFLAGS) $(INCDIR) -c $< -o $@

%.o : %.cpp
	$(CXX) $(CXXFLAGS) $(INCDIR) -c $< -o $@

clean : 
	rm -rf $(patsubst %.cpp,%.o,$(SOURCES)) $(TARGET)
