CC = gcc

LOCAL=${HOME}/.local
CFLAGS = \
	 -pedantic \
	 -Werror \
	 -Wall \
	 -g \
	 -Iinclude \
	 -Wl,-rpath,${LOCAL}/lib \
	 -L${LOCAL}/lib \
	 -I${LOCAL}/include

CFLAGS_OPT = \
	     -O2 \
	     -march=native

LIB = \
      -lm \
      -lgsl \
      -lgslcblas

OBJ = \
      obj/tools.o
MAIN = \
       obj/main.o
LAB2 = \
	obj/lab2.o

ifeq ($(MAKECMDGOALS),test)
OBJ_TEST = $(patsubst obj/%.o,obj_test/%.o,$(OBJ))
-include unit-test/test.mk
else
CFLAGS += $(CFLAGS_OPT)
endif



lab2: obj _lab2

_lab2: $(LAB2) $(OBJ)
	$(CC) $(CFLAGS) $^ -o H2 $(LIB)

obj/%.o: src/%.c
	$(CC) -MMD -c $(CFLAGS) $< -o $@ 

obj:
	mkdir -p obj

clean:
	find -iname "*.o" -exec rm {} \;
	find -iname "*.d" -exec rm {} \;
	rm -f H2 run-test
	rm -rf obj obj_test

.PHONY: clean
