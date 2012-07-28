CC = gcc
CFLAGS = -Wall -O2 -DUSE_DOUBLE -fPIC

CFILES1 = $(shell ls CHAPTER1/*.c)
CFILES2 = $(shell ls CHAPTER2/*.c)
CFILES3 = $(shell ls CHAPTER3/*.c)
CFILES4 = $(shell ls CHAPTER4/*.c)
CFILES5 = $(shell ls CHAPTER5/*.c)
CFILES6 = $(shell ls CHAPTER6/*.c)
CFILES7 = $(shell ls CHAPTER7/*.c)
CFILESU = $(shell ls UTILITY/*.c)

CFILES = $(CFILES1) $(CFILES2) $(CFILES3) $(CFILES4) $(CFILES5) $(CFILES6) \
	 $(CFILES7) $(CFILESU)

OFILES = $(CFILES:.c=.o)

LIBNAME = libnumal.a

default: $(OFILES)
	ar rs $(LIBNAME) $(OFILES)
	$(CC) -shared $(OFILES) -o libnumal.so

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	cd CHAPTER1 && rm -f *.o && cd ..
	cd CHAPTER2 && rm -f *.o && cd ..
	cd CHAPTER3 && rm -f *.o && cd ..
	cd CHAPTER4 && rm -f *.o && cd ..
	cd CHAPTER5 && rm -f *.o && cd ..
	cd CHAPTER6 && rm -f *.o && cd ..
	cd CHAPTER7 && rm -f *.o && cd ..
	cd UTILITY && rm -f *.o && cd ..
	rm -f $(LIBNAME)
	rm -f libnumal.so
