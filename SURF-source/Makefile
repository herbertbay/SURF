# ------------------ Compilation options ------------------------

# Loads math library.
LIBS = -lm
GET = get
# Compiler flags for release
CFLAGS = -Wall -O3 -march=x86-64 -mmmx -msse -msse2 -msse3
# Compiler flags for debugging
# CFLAGS = -Wall -g -pg
# Compiler flags for profiling
# CFLAGS = -Wall -g -O4 -pg
DEPEND= makedepend $(CFLAGS)

CC = g++

# --------------------- Code modules ----------------------------

# Source files
SRCS = main.cpp image.cpp fasthessian.cpp surf.cpp imload.cpp
# Object files
OBJ = main.o image.o fasthessian.o surf.o imload.o
# Object files for shared library
SOBJ = image.o fasthessian.o surf.o
# Definitions
DEFS = image.h fasthessian.h ipoint.h surf.h imload.h

# ------------------------ Rules --------------------------------
#$(SRCS):
#        $(GET) $@

surf.ln: ${OBJ} libSurf.so
	${CC} -o $@ ${CFLAGS} main.o imload.o -L. -lSurf ${LIBS}

libSurf.so: ${SOBJ}
	${CC} ${CFLAGS} -shared -o libSurf.so ${SOBJ} ${LIBS}

libSurf.a: ${SOBJ}
	ar rcs libSurf.a ${SOBJ}

# Implicit rule used by Gnu Make: $(CC) -c $(CPPFLAGS) $(CFLAGS)

imload.o: imload.h imload.cpp image.h
	${CC} -c ${CFLAGS} imload.cpp

surf.o : surf.h surf.cpp
	${CC} -c ${CFLAGS} -fPIC surf.cpp

image.o: image.h image.cpp
	${CC} -c ${CFLAGS} -fPIC image.cpp

fasthessian.o: ${DEFS} fasthessian.cpp
	${CC} -c ${CFLAGS} -fPIC fasthessian.cpp

main.o: ${DEFS} ${SRCS}
	${CC} -c  ${CFLAGS} main.cpp

clean:
	-rm *.o surf.ln libSurf.so*

#depend: $(SRCS)
#        $(DEPEND) $(SRCS)
