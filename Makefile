#CC = /opt/homebrew/bin/gcc-14
CC = gcc
CFLAGS = -O3 -march=native -lm -ffast-math -fopenmp
TARGETS = thrj_000_opt thrj_220_opt

all: $(TARGETS)

thrj_000_opt: thrj_000_opt.c
	$(CC) $(CFLAGS) -o thrj_000_opt thrj_000_opt.c

thrj_220_opt: thrj_220_opt.c
	$(CC) $(CFLAGS) -o thrj_220_opt thrj_220_opt.c

#thrj_array_220: thrj_array_220.c
	#$(CC) $(CFLAGS) -o thrj_array_220 thrj_array_220.c

clean:
	rm -f $(TARGETS)
