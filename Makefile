#
#  (c) 2018-2023 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
#

CFLAGS=-std=c99 -O3 -fno-math-errno -flto -s
WARNINGS=-Wall -Wextra -pedantic -Wshadow -Wpointer-arith -Wcast-qual -Wstrict-prototypes -Wmissing-prototypes -Wconversion -Wredundant-decls -Wmissing-declarations
LIB_STD=-lmpfr
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

all: tsm-std


tsm-%-std: tsm-%.o taylor-ode.o main.o
	$(CC) $(CFLAGS) -o $@ $^ $(LIB_STD)

tsm-std: tsm-bouali-std tsm-burke-shaw-std tsm-genesio-tesi-std tsm-halvorsen-std tsm-isuc-std tsm-lorenz-std tsm-rf-std tsm-rossler-std tsm-rucklidge-std tsm-thomas-std tsm-wimol-banlue-std tsm-yu-wang-std


divergence: divergence.o
	$(CC) $(CFLAGS) -o $@ $< $(LIB_STD)


libad-test: libad-test.o taylor-ode.o
	$(CC) $(CFLAGS) -o $@ $^ $(LIB_STD)

libdual-test: libdual-test.o dual.o
	$(CC) $(CFLAGS) -o $@ $^ $(LIB_STD)

tests: libad-test libdual-test


.PHONY: test clean depclean ctags ctags-system ctags-system-all coverage

test:
	@for x in -.5 0 .5; do \
		./libad-test 20 $$x 1e-15 >/dev/null || exit 1; echo ""; \
		./libdual-test $$x 1e-15 >/dev/null || exit 1; echo ""; \
	done
	@echo "\033[1;37mCore Tests Passed\033[0;37m\n"

test-all: test
	@if ! ./tsm-lorenz-std  6 8 .01 1000  -15.8 -17.48 35.64  10 28 8 3 >/dev/null; then exit 1; fi
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
