CC=mpicc
INC=src/header/
CFLAGS=-O3 -I$(INC)
LIBS=-lm
EXE=reconstruct

COMMON_C=$(wildcard src/*.c)
SERIAL_C=$(wildcard src/serial/*.c)
PARALLEL_C=$(wildcard src/parallel/*.c)

COMMON_O=$(patsubst %.c, %.o, $(COMMON_C))
SERIAL_O=$(patsubst %.c, %.o, $(SERIAL_C))
PARALLEL_O=$(patsubst %.c, %.o, $(PARALLEL_C))

.PHONY: serial
serial: $(SERIAL_O) $(COMMON_O)
	$(CC) $(CFLAGS) $^ -o $(EXE).$@ $(LIBS)

.PHONY: parallel
parallel: $(PARALLEL_O) $(COMMON_O)
	$(CC) $(CFLAGS) $^ -o $(EXE).$@ $(LIBS)

%.o: %.c
	$(CC) $(CFLAGS) -c $^ -o $@ $(LIBS)
.PHONY: clean
clean:
	rm -f $(SERIAL_O) $(PARALLEL_O) $(COMMON_O) $(EXE).* core 
