CC = gcc
CFLAGS = -I${CURDIR}/include -O2
OBJECTS_EU = src/euler.c
OBJECTS_CN = src/crank-nicolson.c

eu: $(OBJECTS_EU)
	$(CC) $(CFLAGS) -o run $(OBJECTS_EU)

cn: $(OBJECTS_CN)
	$(CC) $(CFLAGS) -o run $(OBJECTS_CN)

clean:
	rm -rf run