CC = gcc
CFLAGS = -O3 -march=native -lm -ffast-math -fopenmp
SRC_DIR = src
TARGETS = $(SRC_DIR)/thrj_000_opt $(SRC_DIR)/thrj_220_opt

all: $(TARGETS)

$(SRC_DIR)/thrj_000_opt: $(SRC_DIR)/thrj_000_opt.c
	$(CC) $(CFLAGS) -o $@ $<

$(SRC_DIR)/thrj_220_opt: $(SRC_DIR)/thrj_220_opt.c
	$(CC) $(CFLAGS) -o $@ $<

clean:
	rm -f $(TARGETS)
