CC	  = g++

CXXFLAGS  = -c -Wall -g

OPTIMIZE  = -O2 -m64

LDFLAGS	  = -g

LIBS   	  = #-I/usr/loacl/include -L/usr/local/lib

SOURCES   = readvtklis.cpp

OBJECTS   = $(SOURCES:.cpp=.o)

INCL      = $(SOURCES:.cpp=.h)

EXEC      = readlis

all: $(SOURCES) $(EXEC)

$(EXEC): $(OBJECTS) 
	$(CC) $(OPTIMIZE) $(LDFLAGS) $(OBJECTS)  $(LIBS) -o $@
	rm $(OBJECTS)
.cpp.o:
	$(CC) $(OPTIMIZE) $(CFLAGS) $< -o $@

.PHONY : clean

clean:
	 rm -f $(OBJECTS) $(EXEC)
