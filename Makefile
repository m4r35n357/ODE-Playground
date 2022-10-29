#
#  (c) 2018-2022 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
#

CC=clang -std=c99 -O3 -flto
#CC=gcc -std=c99 -Og -g

CFLAGS=-I. -Wall -Wextra -pedantic -Wshadow -Wpointer-arith -Wcast-qual -Wstrict-prototypes -Wmissing-prototypes -Wconversion -Wredundant-decls -Wmissing-declarations

CFLAGS += -ffp-model=precise
#CFLAGS += -Wunsuffixed-float-constants -frounding-math -fsignaling-nans

LIBS=-lm
GL_LIBS=-lGLEW -lglut -lGLU -lGL

all: ctags tsm tsm-gl h generators h-kerr-gl h-nbody-gl divergence tests

tsm-%.o: tsm-%.c taylor-ode.h real.h
	$(CC) -c -o $@ $< $(CFLAGS)

tsm-%-std: tsm-%.o taylor-ode.o main-tsm.o
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

tsm: tsm-bouali-std tsm-burke-shaw-std tsm-cosmology-std tsm-genesio-tesi-std tsm-halvorsen-std tsm-isuc-std tsm-logistic-std tsm-kom-std tsm-lorenz-std tsm-lotka-volterra-std tsm-nose-hoover-std tsm-rf-std tsm-rossler-std tsm-rucklidge-std tsm-sj-std tsm-sprott-minimal-std tsm-sprott-thomas-std tsm-thomas-std tsm-van-der-pol-std tsm-wimol-banlue-std tsm-yu-wang-std

tsm-%-gl.o: tsm-%.c opengl.h taylor-ode.h real.h
	$(CC) -c -o $@ $< $(CFLAGS)

tsm-%-gl: tsm-%.o taylor-ode.o opengl.o ode-gl.o
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS) $(GL_LIBS)

tsm-gl: tsm-bouali-gl tsm-burke-shaw-gl tsm-cosmology-gl tsm-genesio-tesi-gl tsm-halvorsen-gl tsm-isuc-gl tsm-logistic-gl tsm-kom-gl tsm-lorenz-gl tsm-lotka-volterra-gl tsm-nose-hoover-gl tsm-rf-gl tsm-rossler-gl tsm-rucklidge-gl tsm-sj-gl tsm-sprott-minimal-gl tsm-sprott-thomas-gl tsm-thomas-gl tsm-van-der-pol-gl tsm-wimol-banlue-gl tsm-yu-wang-gl

h-%.o: h-%.c symplectic.h dual.h real.h
	$(CC) -c -o $@ $< $(CFLAGS)

h-%-std: h-%.o symplectic.o dual.o
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

h: h-analysis-std h-newton-std h-spring-std

h-kerr-gen-light: h-kerr-gen-light.c real.h
	$(CC) -o $@ $< $(CFLAGS) $(LIBS)

h-kerr-gen-particle.o: h-kerr-gen-particle.c h-kerr.h dual.h real.h
	$(CC) -c -o $@ $< $(CFLAGS)

h-kerr-gen-particle: h-kerr-gen-particle.o h-kerr.o dual.o
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

generators: h-kerr-gen-light h-kerr-gen-particle

h-kerr-gl.o: h-kerr-gl.c opengl.h h-kerr.h symplectic.h dual.h real.h
	$(CC) -c -o $@ $< $(CFLAGS)

h-kerr-gl: h-kerr-gl.o opengl.o symplectic.o dual.o h-kerr.o
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS) $(GL_LIBS)

h-nbody-gl.o: h-nbody-gl.c opengl.h h-nbody.h symplectic.h real.h
	$(CC) -c -o $@ $< $(CFLAGS)

h-nbody-gl: h-nbody-gl.o opengl.o symplectic.o h-nbody.o
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS) $(GL_LIBS)

divergence: divergence.c real.h
	$(CC) -o $@ $< $(CFLAGS) $(LIBS)

libad-test.o: libad-test.c taylor-ode.h ad.h real.h
	$(CC) -c -o $@ $< $(CFLAGS)

libad-test: libad-test.o taylor-ode.o ad.o
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

libdual-test.o: libdual-test.c dual.h real.h
	$(CC) -c -o $@ $< $(CFLAGS)

libdual-test: libdual-test.o dual.o
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

tests: libad-test libdual-test

kerr-image: kerr-image.c
	$(CC) -o $@ $< $(CFLAGS) $(LIBS)

.PHONY: clean ctags

ctags:
	/usr/bin/ctags *.h *.c

clean:
	rm -f *.o *~ core
	rm -f *-std
	rm -f *-gl
	rm -f divergence
	rm -f h-kerr-gen-light h-kerr-gen-particle
	rm -f libad-test libdual-test
