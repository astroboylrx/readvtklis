CC	  = g++

CFLAGS	  = -c -Wall -g

OPTIMIZE  = -O2

LDFLAGS	  = -g

LIBS   	  = #-I/usr/loacl/include -L/usr/local/lib

SOURCES   = readvtklis.cpp

OBJECTS   = $(SOURCES:.cpp=.o)

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
