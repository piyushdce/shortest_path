# shortest_path
Dijkstra’s algorithm is widely used to calculate single source shortest path in graphs having non negative weights. In this project Dijkstra’s algorithm is implemented in C programming using efficient data structures and is tested on the US continental map having 87k vertices and 121k edges. The algorithm takes only 2.2 ms on average to return the shortest distance and the path between two nodes in the graph.  
An improved version of Dijkstra’s algorithm (A star) directs its search towards the destination based on an optimistic heuristic distance. The algorithm when implemented improves the shortest distance calculation time by 2.5x by reducing the number of extract min and decrease key operations by 5x.
Usage: ./SP_Dijkstra.c <Map File> <query File>
  Where Map File is usa.txt or your own Map file which contain the node number followed by the two vertices it connects and query file is the list of query pairs (node1 node2). usa100.txt has 100 query items.
