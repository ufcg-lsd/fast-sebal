CC=g++
CFLAGS = -std=c++14
CLTIFF = -ltiff

DEPS = src/*.h

SRCS = $(shell find -name *.cpp)
OBJS := $(addsuffix .o,$(basename $(SRCS)))

all: run

src/%.o: src/%.cpp $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS) $(CLTIFF)
	
run: $(SRCS)
	$(CC) -o $@ $^ $(CFLAGS) $(CLTIFF)

clean: 
	rm -f run
	rm -rf *.o
	rm -rf src/*.o