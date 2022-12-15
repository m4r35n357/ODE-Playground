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
  CC=gcc -std=c99 -Og -g -pg --coverage
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

all: tsm tsm-gl h generators h-kerr-std h-kerr-gl h-nbody-std h-nbody-gl divergence tests


tsm-%-std: tsm-%.o taylor-ode.o main-tsm.o
	$(CC) -o $@ $^ $(LIB_M) $(STRIP)

tsm: tsm-bouali-std tsm-burke-shaw-std tsm-cosmology-std tsm-genesio-tesi-std tsm-halvorsen-std tsm-isuc-std tsm-logistic-std tsm-kom-std tsm-lorenz-std tsm-lotka-volterra-std tsm-nose-hoover-std tsm-rf-std tsm-rossler-std tsm-rucklidge-std tsm-sprott-minimal-std tsm-sprott-thomas-std tsm-thomas-std tsm-van-der-pol-std tsm-wimol-banlue-std tsm-yu-wang-std


tsm-%-gl: tsm-%.o taylor-ode.o opengl.o ode-gl.o
	$(CC) -o $@ $^ $(LIB_M) $(LIB_GL) $(STRIP)

tsm-gl: tsm-bouali-gl tsm-burke-shaw-gl tsm-cosmology-gl tsm-genesio-tesi-gl tsm-halvorsen-gl tsm-isuc-gl tsm-logistic-gl tsm-kom-gl tsm-lorenz-gl tsm-lotka-volterra-gl tsm-nose-hoover-gl tsm-rf-gl tsm-rossler-gl tsm-rucklidge-gl tsm-sprott-minimal-gl tsm-sprott-thomas-gl tsm-thomas-gl tsm-van-der-pol-gl tsm-wimol-banlue-gl tsm-yu-wang-gl


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


.PHONY: test clean depclean ctags ctags-system ctags-system-all coverage

test:
	@for x in -.5 0 .5; do \
		./libad-test 20 $$x 1e-15 >/dev/null || exit 1; ./libdual-test $$x 1e-15 >/dev/null || exit 1; \
	done
	@if ! ./tsm-lorenz-std  6 8 .01 1000  -15.8 -17.48 35.64  10 28 8 3 >/dev/null; then exit 1; fi
	@if ! ./h-newton-std  6 10 1 1000  1 12 .6 >/dev/null; then exit 1; fi
	@if ! ./h-kerr-std 6 8 0.010 1000 _ 0.8 1.0 0.9455050956749083 1.434374509531738 1.0 7.978759958927879 12.0 63.0 >/dev/null; then exit 1; fi
	@if ! ./h-nbody-std 6 6 0.010 1000 0.05 100.0 0.0 0.0 0.0 0.0 0.0 0.0 2.0 0.0 4.5 0.4 -0.2 0.0 1.8 3.0 -6.0 0.0 -0.4 0.0 -2.0 1.0 5.0 3.0 0.0 -0.2 0.0 5.8 -0.2 4.0 0.0 -4.0 0.1 -3.6 0.0 0.2 3.0 -4.0 0.0 -0.1 0.0 -0.2 -2.6 3.0 8.0 0.0 -0.3 0.0 2.0 -0.2 4.0 0.0 4.0 -0.2 -4.8 0.0 -0.2 >/dev/null; then exit 1; fi
	@if ! ./h-kerr-gen-particle 1.0e-9 4.0 12.0 63.0 0.8 >/dev/null; then exit 1; fi
	@if ! ./h-kerr-gen-light 3.0 0.8 >/dev/null; then exit 1; fi

ctags:
	@/usr/bin/ctags *.h *.c

ctags-system:
	@/usr/bin/ctags --c-kinds=+p --fields=+iaS --extras=+q /usr/include/*.h /usr/include/GL/*.h *.h *.c

ctags-system-all:
	@/usr/bin/ctags -R --c-kinds=+p --fields=+iaS --extras=+q /usr/include .

clean:
	@rm -f *.o *.gcda *.gcno *~ core *-std *-gl h-kerr-gen-light h-kerr-gen-particle divergence libad-test libdual-test

depclean: clean
	@rm -f *.d

C_OUT_DIR=coverage-out
C_OUT_FILE=coverage.info
coverage: test
	@lcov --capture --directory . --output-file $(C_OUT_FILE)
	@genhtml $(C_OUT_FILE) --output-directory $(C_OUT_DIR)
	@echo "C:      file://$(CURDIR)/$(C_OUT_DIR)/index.html"

-include *.d
