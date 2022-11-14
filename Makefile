#
#  (c) 2018-2022 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
#

CFLAGS=-Wall -Wextra -pedantic -Wshadow -Wpointer-arith -Wcast-qual -Wstrict-prototypes -Wmissing-prototypes -Wconversion -Wredundant-decls -Wmissing-declarations
LIB_M=-lm
LIB_GL=-lGLEW -lglut -lGLU -lGL
STRIP=-s

ifeq ($(CCC),gcc)
  CC=gcc -std=c99 -O3 -flto
  CFLAGS += -Wunsuffixed-float-constants -frounding-math -fsignaling-nans
else ifeq ($(CCC),dbg)
  CC=gcc -std=c99 -Og -g
  CFLAGS += -Wunsuffixed-float-constants -frounding-math -fsignaling-nans
  STRIP=
else ifeq ($(CCC),clang)
  CC=clang -std=c99 -O3 -flto
  CFLAGS += -ffp-model=precise
else
  CC=clang -std=c99 -O3 -flto
  CFLAGS += -ffp-model=precise
endif

%.o: %.c
	$(CC) -MT $@ -MMD -MP -c -o $@ $< $(CFLAGS)

all: ctags tsm tsm-gl h generators h-kerr-std h-kerr-gl h-nbody-std h-nbody-gl divergence tests


tsm-%-std: tsm-%.o taylor-ode.o main-tsm.o
	$(CC) -o $@ $^ $(LIB_M) $(STRIP)

tsm: tsm-bouali-std tsm-burke-shaw-std tsm-cosmology-std tsm-genesio-tesi-std tsm-halvorsen-std tsm-isuc-std tsm-logistic-std tsm-kom-std tsm-lorenz-std tsm-lotka-volterra-std tsm-nose-hoover-std tsm-rf-std tsm-rossler-std tsm-rucklidge-std tsm-sj-std tsm-sprott-minimal-std tsm-sprott-thomas-std tsm-thomas-std tsm-van-der-pol-std tsm-wimol-banlue-std tsm-yu-wang-std


tsm-%-gl: tsm-%.o taylor-ode.o opengl.o ode-gl.o
	$(CC) -o $@ $^ $(LIB_M) $(LIB_GL) $(STRIP)

tsm-gl: tsm-bouali-gl tsm-burke-shaw-gl tsm-cosmology-gl tsm-genesio-tesi-gl tsm-halvorsen-gl tsm-isuc-gl tsm-logistic-gl tsm-kom-gl tsm-lorenz-gl tsm-lotka-volterra-gl tsm-nose-hoover-gl tsm-rf-gl tsm-rossler-gl tsm-rucklidge-gl tsm-sj-gl tsm-sprott-minimal-gl tsm-sprott-thomas-gl tsm-thomas-gl tsm-van-der-pol-gl tsm-wimol-banlue-gl tsm-yu-wang-gl


h-%-std: h-%.o symplectic.o dual.o
	$(CC) -o $@ $^ $(LIB_M) $(STRIP)

h: h-analysis-std h-newton-std h-spring-std


h-kerr-gen-light: h-kerr-gen-light.o
	$(CC) -o $@ $< $(LIB_M) $(STRIP)

h-kerr-gen-particle: h-kerr-gen-particle.o h-kerr.o dual.o
	$(CC) -o $@ $^ $(LIB_M) $(STRIP)

generators: h-kerr-gen-light h-kerr-gen-particle


h-kerr-std: symplectic.o dual.o h-kerr.o main-kerr.o
	$(CC) -o $@ $^ $(LIB_M) $(STRIP)

h-kerr-gl: symplectic.o dual.o h-kerr.o opengl.o h-kerr-gl.o
	$(CC) -o $@ $^ $(LIB_M) $(LIB_GL) $(STRIP)


h-nbody-std: symplectic.o h-nbody.o main-nbody.o
	$(CC) -o $@ $^ $(LIB_M) $(LIB_GL) $(STRIP)

h-nbody-gl: symplectic.o h-nbody.o opengl.o h-nbody-gl.o
	$(CC) -o $@ $^ $(LIB_M) $(LIB_GL) $(STRIP)


divergence: divergence.o
	$(CC) -o $@ $< $(LIB_M) $(STRIP)


libad-test: libad-test.o taylor-ode.o ad.o
	$(CC) -o $@ $^ $(LIB_M) $(STRIP)

libdual-test: libdual-test.o dual.o
	$(CC) -o $@ $^ $(LIB_M) $(STRIP)

tests: libad-test libdual-test


kerr-image: kerr-image.o
	$(CC) -o $@ $< $(LIB_M) $(STRIP)


.PHONY: clean depclean ctags

ctags:
	@/usr/bin/ctags *.h *.c

clean:
	@rm -f *.o *~ core *-std *-gl h-kerr-gen-light h-kerr-gen-particle divergence libad-test libdual-test

depclean: clean
	@rm -f *.d

-include *.d
