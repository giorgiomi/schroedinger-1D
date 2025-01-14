CC = gcc
CFLAGS = -I${CURDIR}/include -O2
OBJECTS_EU = src/euler.c src/functions.c
OBJECTS_CN = src/crank-nicolson.c src/functions.c
OBJECTS_TS = src/trotter-suzuki.c src/functions.c

eu: $(OBJECTS_EU)
	$(CC) $(CFLAGS) -o eu.x $(OBJECTS_EU)

cn: $(OBJECTS_CN)
	$(CC) $(CFLAGS) -o cn.x $(OBJECTS_CN)

ts: $(OBJECTS_TS)
	$(CC) $(CFLAGS) -o ts.x $(OBJECTS_TS)

clean:
	rm -rf run