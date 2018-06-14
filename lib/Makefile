CC := gcc
CP := g++
LIBFLAGS := -L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_intel_ilp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl
INCLUDEFLAGS := -DMKL_ILP64 -m64 -I${MKLROOT}/include 
SHAREFLAGS := -O3 -shared -fPIC

all: reverse_regression.so linear_regression.so
linreg: linear_regression.so
revreg: reverse_regression.so

linear_regression.so: linear_regression.c
	$(CC) $< $(SHAREFLAGS) $(INCLUDEFLAGS) -o $@

reverse_regression.so: reverse_regression.c
	$(CC) $< $(LIBFLAGS) $(SHAREFLAGS) $(INCLUDEFLAGS) -o $@

clean:
	rm -rf linear_regression.so reverse_regression.so 
