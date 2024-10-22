TARGET = ref.exe

OBJECT = main.o func.o

CC = g++

CFLAGS = -c

all: $(OBJECT)
	$(CC) $(OBJECT) -o $(TARGET)	
main.o: main.cpp
	$(CC) $(CFLAGS) main.cpp
func.o: func.cpp
	$(CC) $(CFLAGS) func.cpp
clean:
	rm -f $(OBJECT) $(TARGET)
