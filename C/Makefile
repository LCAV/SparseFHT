CC=gcc
CDEBUG=-g
CFLAGS=-I. -Wall 
DEPS = gf2.h fht.h sfht.h common.h
OBJS = gf2.o fht.o sfht.o common.o

UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Linux)
	CFLAGS += -lm
endif

%.o: %.c $(DEPS)
	    $(CC) $(CDEBUG) -c -o $@ $< $(CFLAGS)
			    
all: gf2_test fht_test sfht_test sfht_det_test common_test benchmark SFHT_ErrorSim ThresholdBehavior

benchmark: $(OBJS) benchmark.o
	$(CC) $(CFLAGS) -o benchmark benchmark.o $(OBJS) -I. -lm -O3

common_test: $(OBJS) common_test.o
	$(CC) $(CDEBUG) $(CFLAGS) -o common_test common_test.o $(OBJS) -I. 

sfht_test: $(OBJS) sfht_test.o
	$(CC) $(CDEBUG) $(CFLAGS) -o sfht_test sfht_test.o $(OBJS) -I. -lm

sfht_det_test: $(OBJS) sfht_det_test.o
	$(CC) $(CDEBUG) $(CFLAGS) -o sfht_det_test sfht_det_test.o $(OBJS) -I. -lm

SFHT_ErrorSim: $(OBJS) SFHT_ErrorSim.o
	$(CC) $(CDEBUG) $(CFLAGS) -o SFHT_ErrorSim SFHT_ErrorSim.o $(OBJS) -I. -lm

ThresholdBehavior: $(OBJS) ThresholdBehavior.o
	$(CC) $(CDEBUG) $(CFLAGS) -o ThresholdBehavior ThresholdBehavior.o $(OBJS) -I. -lm
	      
fht_test: $(OBJS) fht_test.o
	$(CC) $(CDEBUG) $(CFLAGS) -o fht_test fht_test.o $(OBJS) -I.
	      
gf2_test: $(OBJS) gf2_test.o
	$(CC) $(CDEBUG) $(CFLAGS) -o gf2_test gf2_test.o $(OBJS) -I.

clean:
	  rm -rf *.o

