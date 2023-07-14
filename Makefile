#
#  (c) 2018-2023 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
#

CFLAGS=-std=c99 -O3 -fno-math-errno -flto -s
WARNINGS=-Wall -Wextra -pedantic -Wshadow -Wpointer-arith -Wcast-qual -Wstrict-prototypes -Wmissing-prototypes -Wconversion -Wredundant-decls -Wmissing-declarations
LIB_STD=-lm
LIB_GL=-lGLEW -lglut -lGLU -lGL

ifeq ($(CCC),gcc)  # fast option
  CC=/usr/bin/gcc
  WARNINGS += -Wunsuffixed-float-constants
else ifeq ($(CCC),cov)  # coverage
  CC=/usr/bin/gcc
  CFLAGS=-std=c99 -O0 -g --coverage
  WARNINGS += -Wunsuffixed-float-constants
else ifeq ($(CCC),prof)  # profiling with GCC
  CC=/usr/bin/gcc
  CFLAGS=-std=c99 -O0 -g -pg
  WARNINGS += -Wunsuffixed-float-constants
else ifeq ($(CCC),gpt)  # profiling with Clang and Google tools
  CC=/usr/bin/clang
  CFLAGS=-std=c99 -O0 -g
  LIB_STD += -lprofiler
else ifeq ($(CCC),clang)  # fast option
  CC=/usr/bin/clang
else  # default for IDEs and git commits
  CC=/usr/bin/clang
  CFLAGS=-std=c99 -O0 -g
endif

%.o: %.c
	$(CC) $(CFLAGS) -MT $@ -MMD -MP -c -o $@ $< $(WARNINGS)

all: tsm-std tsm-gl hamiltonian generators h-kerr-std h-kerr-gl h-nbody-std h-nbody-gl divergence tests


tsm-%-std: tsm-%.o taylor-ode.o main-tsm.o
	$(CC) $(CFLAGS) -o $@ $^ $(LIB_STD)

tsm-std: tsm-bouali-std tsm-burke-shaw-std tsm-genesio-tesi-std tsm-halvorsen-std tsm-isuc-std tsm-lorenz-std tsm-rf-std tsm-rossler-std tsm-rucklidge-std tsm-thomas-std tsm-wimol-banlue-std tsm-yu-wang-std


tsm-%-gl: tsm-%.o taylor-ode.o opengl.o ode-gl.o
	$(CC) $(CFLAGS) -o $@ $^ $(LIB_STD) $(LIB_GL)

tsm-gl: tsm-bouali-gl tsm-burke-shaw-gl tsm-genesio-tesi-gl tsm-halvorsen-gl tsm-isuc-gl tsm-lorenz-gl tsm-rf-gl tsm-rossler-gl tsm-rucklidge-gl tsm-thomas-gl tsm-wimol-banlue-gl tsm-yu-wang-gl


h-%-std: h-%.o symplectic.o dual.o
	$(CC) $(CFLAGS) -o $@ $^ $(LIB_STD)

hamiltonian: h-analysis-std h-newton-std h-spring-std


h-kerr-gen-light: h-kerr-gen-light.o
	$(CC) $(CFLAGS) -o $@ $< $(LIB_STD)

h-kerr-gen-particle: h-kerr-gen-particle.o h-kerr.o dual.o
	$(CC) $(CFLAGS) -o $@ $^ $(LIB_STD)

generators: h-kerr-gen-light h-kerr-gen-particle


h-kerr-std: symplectic.o dual.o h-kerr.o main-kerr.o
	$(CC) $(CFLAGS) -o $@ $^ $(LIB_STD)

h-kerr-gl: symplectic.o dual.o h-kerr.o opengl.o h-kerr-gl.o
	$(CC) $(CFLAGS) -o $@ $^ $(LIB_STD) $(LIB_GL)


h-nbody-std: symplectic.o h-nbody.o main-nbody.o
	$(CC) $(CFLAGS) -o $@ $^ $(LIB_STD) $(LIB_GL)

h-nbody-gl: symplectic.o h-nbody.o opengl.o h-nbody-gl.o
	$(CC) $(CFLAGS) -o $@ $^ $(LIB_STD) $(LIB_GL)


divergence: divergence.o
	$(CC) $(CFLAGS) -o $@ $< $(LIB_STD)


libad-test: libad-test.o taylor-ode.o
	$(CC) $(CFLAGS) -o $@ $^ $(LIB_STD)

libdual-test: libdual-test.o dual.o
	$(CC) $(CFLAGS) -o $@ $^ $(LIB_STD)

tests: libad-test libdual-test


kerr-image: kerr-image.o
	$(CC) $(CFLAGS) -o $@ $< $(LIB_STD)


.PHONY: test clean depclean ctags ctags-system ctags-system-all coverage

test:
	@for x in -.5 0 .5; do \
		./libad-test 20 $$x 1e-15 >/dev/null || exit 1; echo ""; \
		./libdual-test $$x 1e-15 >/dev/null || exit 1; echo ""; \
	done
	@echo "\033[1;37mCore Tests Passed\033[0;37m\n"

test-all: test
	@if ! ./tsm-lorenz-std  6 8 .01 1000  -15.8 -17.48 35.64  10 28 8 3 >/dev/null; then exit 1; fi
	@if ! ./h-analysis-std 6 4 1.0 1 >/dev/null; then exit 1; fi
	@if ! ./h-newton-std  6 10 1 1000  1 12 .6 >/dev/null; then exit 1; fi
	@if ! ./h-kerr-std 6 8 0.010 1000 0.8 1.0 0.9455050956749083 1.434374509531738 1.0 7.978759958927879 12.0 63.0 >/dev/null; then exit 1; fi
	@if ! ./h-nbody-std 6 6 0.010 1000 0.05 100.0 0.0 0.0 0.0 0.0 0.0 0.0 2.0 0.0 4.5 0.4 -0.2 0.0 1.8 3.0 -6.0 0.0 -0.4 0.0 -2.0 1.0 5.0 3.0 0.0 -0.2 0.0 5.8 -0.2 4.0 0.0 -4.0 0.1 -3.6 0.0 0.2 3.0 -4.0 0.0 -0.1 0.0 -0.2 -2.6 3.0 8.0 0.0 -0.3 0.0 2.0 -0.2 4.0 0.0 4.0 -0.2 -4.8 0.0 -0.2 >/dev/null; then exit 1; fi
	@if ! ./h-kerr-gen-particle 1.0e-9 4.0 12.0 63.0 0.8 >/dev/null; then exit 1; fi
	@if ! ./h-kerr-gen-light 3.0 0.8 >/dev/null; then exit 1; fi
	@echo "\n\033[1;37mSanity Tests Passed\033[0;37m"

ctags:
	@/usr/bin/ctags *.h *.c

ctags-system:
	@/usr/bin/ctags --c-kinds=+p --fields=+iaS --extras=+q /usr/include/*.h /usr/include/GL/*.h *.h *.c

ctags-system-all:
	@/usr/bin/ctags -R --c-kinds=+p --fields=+iaS --extras=+q /usr/include .

clean:
	@rm -rf *.o *.gcda *.gcno *-std *-gl h-kerr-gen-light h-kerr-gen-particle divergence libad-test libdual-test \
		coverage* gmon.out

depclean: clean
	@rm -f *.d

C_OUT_DIR=coverage-out
C_OUT_FILE=coverage.info
coverage: test
	@lcov --capture --directory . --output-file $(C_OUT_FILE)
	@genhtml $(C_OUT_FILE) --output-directory $(C_OUT_DIR)
	@echo "C:      file://$(CURDIR)/$(C_OUT_DIR)/index.html"

-include *.d
