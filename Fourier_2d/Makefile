TARGET_MAIN = ref.exe
TARGET_CONV = conv.exe

OBJECTS_MAIN = main.o func.o
OBJECTS_CONV = convergence.o func.o

CC = g++

CFLAGS = -c

all: $(OBJECTS_MAIN)
	$(CC) $(OBJECTS_MAIN) -o $(TARGET_MAIN)	
conv: $(OBJECTS_CONV)
	$(CC) $(OBJECTS_CONV) -o $(TARGET_CONV)
main.o: main.cpp
	$(CC) $(CFLAGS) main.cpp
func.o: func.cpp
	$(CC) $(CFLAGS) func.cpp
convergence.o: convergence.cpp
	$(CC) $(CFLAGS) convergence.cpp
clean:
	rm -f $(OBJECTS_MAIN) $(TARGET_MAIN)
	rm -f $(OBJECTS_CONV) $(TARGET_CONV)

