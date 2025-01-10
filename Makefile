CC = gcc
CFLAGS = -I${CURDIR}/include -O2
OBJECTS = src/main.c

run: $(OBJECTS)
	$(CC) $(CFLAGS) -o run $(OBJECTS)
clean:
	rm -rf run