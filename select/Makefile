BIN=select.bin
SRC=select.c

all : ${BIN} s.bin best_sphere.bin

${BIN} : ${SRC}
	gcc $^ -O2 -lm -o $@

s.bin : s.c
	gcc -O2 -lm s.c -o s.bin

best_sphere.bin : best_sphere.c
	gcc -O2 -lm best_sphere.c -o best_sphere.bin
