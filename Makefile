.PHONY: build clean

build: filtru

filtru: main.o functions.o
	mpicc $^ -o $@

main.o: main.c
	mpicc -c $<

functions.o: functions.c
	mpicc -c $<

clean:
	rm -f *.o filtru
