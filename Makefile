CC=gcc
DEPS= average.h 
OBJ= main.o average.o
CFLAGS=-std=gnu99

%.o: %.c $(DEPS)
	$(CC) $(CFLAGS) -c -o $@ $<

main: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS)
