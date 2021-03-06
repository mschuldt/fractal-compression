
OPT = -fopenmp -msse4 -O2 counters.cpp -std=c++11

OBJ = IFSTransform.o\
	Image.o\
	Decoder.o\
	Encoder.o\
	QuadTreeEncoder.o\
        count_ops.o


all: fractal

IFSTransform.o: IFSTransform.h IFSTransform.cpp
	g++ $(OPT) -c IFSTransform.cpp

Image.o: Image.h Image.cpp
	g++ $(OPT) -c Image.cpp

Decoder.o: Decoder.h Decoder.cpp
	g++ $(OPT) -c Decoder.cpp

Encoder.o: Encoder.h Encoder.cpp
	g++ $(OPT) -c Encoder.cpp

count_ops.o: count_ops.cpp
	g++ $(OPT) -c count_ops.cpp

QuadTreeEncoder.o: QuadTreeEncoder.h QuadTreeEncoder.cpp
	g++ $(OPT) -c QuadTreeEncoder.cpp

fractal: $(OBJ) main.cpp
	g++ $(OPT) -o fractal $(OBJ) main.cpp

check:
	cd test; ./test.sh

clean:
	rm *.o
	rm fractal*


