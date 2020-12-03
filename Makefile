INCLUDE=../dSFMT-src-2.2.3
# si cambiaron los componentes recompila, si no, no hace falta.
# no hace falta borrar el ejecutable

Pcorona: Xcorona.h Pcorona.c
	gcc -O3 -o Pcorona -lm -IINCLUDE -LINCLUDE -DDSFMT_MEXP=19937  Pcorona.c ${INCLUDE}/dSFMT.o

Icorona: Xcorona.h Icorona.c
	gcc -O3 -o Icorona -lm -IINCLUDE -LINCLUDE -DDSFMT_MEXP=19937  Icorona.c ${INCLUDE}/dSFMT.o

Acorona: Xcorona.h Acorona.c
	gcc -O3 -o Acorona -lm -IINCLUDE -LINCLUDE -DDSFMT_MEXP=19937  Acorona.c ${INCLUDE}/dSFMT.o

corona: corona.h corona.c
	gcc -O3 -o corona -lm -IINCLUDE -LINCLUDE -DDSFMT_MEXP=19937  corona.c ${INCLUDE}/dSFMT.o
