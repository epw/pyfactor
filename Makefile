CFLAGS = -g -Wall `pkg-config --cflags python` --shared -fpic
LIBS = -lm `pkg-config --libs python`

all: pyfactor.so

PYFACTOR_OBJS = pyfactor.o error.o quote.o
pyfactor.so: $(PYFACTOR_OBJS)
	$(CC) $(CFLAGS) -o pyfactor.so $(PYFACTOR_OBJS) $(LIBS)

clean:
	rm -f *.o *.so
