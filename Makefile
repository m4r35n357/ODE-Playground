
CC=gcc -std=c99
CFLAGS=-I. -Wall -Wextra -pedantic

LIBS=-lm
GL_LIBS=-lGLEW -lglut -lGLU -lGL

TSM_DEPS = taylor-ode.h real.h
TSM_OBJ = taylor-ode.o main-tsm.o

tsm-%.o: tsm-%.c $(TSM_DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

tsm-%-std: tsm-%.o $(TSM_OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

tsm: tsm-halvorsen-std tsm-lorenz-std tsm-thomas-std tsm-rf-std tsm-rossler-std

TSM_GL_DEPS = opengl.h taylor-ode.h real.h
TSM_GL_OBJ = taylor-ode.o opengl.o ode-gl.o

tsm-%-gl.o: tsm-%.c $(TSM_GL_DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

tsm-%-gl: tsm-%.o $(TSM_GL_OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS) $(GL_LIBS)

tsm-gl: tsm-halvorsen-gl tsm-lorenz-gl tsm-thomas-gl tsm-rf-gl tsm-rossler-gl

H_DEPS = symplectic.h real.h
H_OBJ = symplectic.o dual.o

h-%.o: h-%.c $(H_DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

h-%-std: h-%.o $(H_OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

h: h-analysis-std h-newton-std h-spring-std

h-kerr-gen-light: h-kerr-gen-light.c real.h
	$(CC) -o $@ $< $(CFLAGS) $(LIBS)

h-kerr-gen-particle.o: h-kerr-gen-particle.c real.h
	$(CC) -c -o $@ $< $(CFLAGS)

h-kerr-gen-particle: h-kerr-gen-particle.o dual.o h-kerr.o
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

H_KERR_GL_DEPS = opengl.h h-kerr.h symplectic.h dual.h real.h
H_KERR_GL_OBJ = opengl.o symplectic.o dual.o h-kerr.o

h-kerr-gl.o: h-kerr-gl.c $(H_KERR_GL_DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

h-kerr-gl: h-kerr-gl.o $(H_KERR_GL_OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS) $(GL_LIBS)

H_NBODY_GL_DEPS = opengl.h h-nbody.h symplectic.h real.h
H_NBODY_GL_OBJ = opengl.o symplectic.o h-nbody.o

h-nbody-gl.o: h-nbody-gl.c $(H_NBODY_GL_DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

h-nbody-gl: h-nbody-gl.o $(H_NBODY_GL_OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS) $(GL_LIBS)

divergence: divergence.c real.h
	$(CC) -o $@ $< $(CFLAGS) $(LIBS)

AD_TEST_DEPS = taylor-ode.h ad.h real.h
AD_TEST_OBJ = taylor-ode.o ad.o

libad-test.o: libad-test.c $(AD_TEST_DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

libad-test: libad-test.o $(AD_TEST_OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

DUAL_TEST_DEPS = dual.h real.h
DUAL_TEST_OBJ = dual.o

libdual-test.o: libdual-test.c $(DUAL_TEST_DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

libdual-test: libdual-test.o $(DUAL_TEST_OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

all: tsm tsm-gl h h-kerr-gen-light h-kerr-gen-particle h-kerr-gl h-nbody-gl divergence libad-test libdual-test

default: all

.PHONY: clean

clean:
	rm -f *.o *~ core
	rm -f *-std
	rm -f *-gl
	rm -f divergence
	rm -f h-kerr-gen-light h-kerr-gen-particle
	rm -f libad-test libdual-test
