CXX	  = g++

OPTIMIZE  = -O2 -m64 -Wall -g

SOURCES   = readvtklis.cpp readvtk.cpp readlis.cpp fop.cpp
OBJECTS   = $(SOURCES:.cpp=.o)
INCL      = readvtk.h readlis.h fop.h

CXXFLAGS  = #-I/opt/local/include/gcc48
LDFLAGS   = #-L/opt/local/lib/gcc48
LIBS      = -lm

EXEC      = readvtklis

all: $(SOURCES) $(EXEC)

$(EXEC): $(OBJECTS) 
	$(CXX) $(OPTIMIZE) $(LDFLAGS) $(OBJECTS) $(LIBS) -o $@
	rm $(OBJECTS)
.cpp.o:
	$(CXX) $(OPTIMIZE) $(CXXFLAGS) -c $< -o $@

.PHONY : clean

clean:
	 rm -f $(OBJECTS) $(EXEC)
