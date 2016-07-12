CC=mpicc

all: serial pagerank

serial: pagerank_serial.c
				$(CC) pagerank_serial.c -o serial

pagerank: pagerank.c
				$(CC) pagerank.c -o pagerank

clean:
				rm pagerank serial