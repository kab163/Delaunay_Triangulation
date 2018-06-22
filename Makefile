SRC = DTcode.C logging.C
OBJ=$(SRC:.C=.o)

all: $(OBJ)
	g++ -w -g -I. $(OBJ) -o delt

.C.o: $<
	g++ -w -g -I. -c $<

run:
	./delt

clean:
	rm *.o delt
