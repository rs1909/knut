noinst_LIBRARIES = liblaarpack.a
liblaarpack_a_SOURCES = laarpack.c dlamch.c cblaswrap.h f2c.h

INCLUDES = -I../ATLAS/include

dlamch.o: dlamch.c
	$(CC) -c @NOOPT_CFLAGS@ $(INCLUDES) dlamch.c
