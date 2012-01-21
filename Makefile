CFLAGS = -g -Wall `pkg-config --cflags python` --shared -fpic
LIBS = -lm -lgmp `pkg-config --libs python`

all: factor.so

FACTOR_OBJS = factor.o
factor.so: $(FACTOR_OBJS)
	$(CC) $(CFLAGS) -o factor.so $(FACTOR_OBJS) $(LIBS)

clean:
	rm -f *.o *.so
