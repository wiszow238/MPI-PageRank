Typing "make" will compile both the parallel and serial program. 

Run commands:

To run the serial program type:
$ ./serial <graphFile>

To run the parallel program type:
$ mpirun -machinefile machines -np <# of Partitions> pagerank <graph file> <Partition File>

Example:
$ mpirun -machinefile machines -np 4 pagerank 4M-graph.txt 4M-graph.txt.part.4

The results will be printed out in "pagerank.result".

I start my time right before I start to loop and calculate the normal values. The start time can be found on line 463. Before line 463, the data structures are created and populated. Once everthing is populated then the computations can be performed, which is where the timer is started. 

The timer ends right when the change in normal values is less than 10^{-5}. The finish can be found on line 572.