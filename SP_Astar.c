#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <stdbool.h>

/******************************
*[Piyush Agarwal]
*[SP_Astar.c]: [7/26/18].
*[agarw145@purdue.edu]
*Compiling: [gcc -Werror -lm -Wall -O3 SP_Astar.c -o Astar]
*******************************/

struct heapnode
	{
	int value;
	int vertex;
	};

struct MinHeap
	{
    int size;
    struct heapnode *array;
	};

struct AdjListNode 		// Adj list node
	{
    int adjNode;
	int weight;
    struct AdjListNode* next;
	};

struct AdjList 			// adjacency list source
	{
	int parent;
	int distance;
	int node_heap;
	short int visited;
    struct AdjListNode *head; 
	};

struct Graph 			// structure for graph
	{
    int V;
    struct AdjList* array;
	struct AdjListNode* edge_array;
	};

// Create graph of V vertices.
struct Graph* createGraph(int V, int E)
	{
    struct Graph* graph = (struct Graph*) malloc(sizeof(struct Graph));
    graph->V = V;
    graph->array = (struct AdjList*) malloc(V * sizeof(struct AdjList));
    graph->edge_array = (struct AdjListNode*) malloc(2*E*sizeof(struct AdjListNode));
 
    int i;
    for (i = 0; i < V; i++)
        {
		graph->array[i].parent = INT_MAX;
		graph->array[i].distance = INT_MAX;
		graph->array[i].node_heap = -1;
		graph->array[i].visited = 0;
		graph->array[i].head = NULL;
		}
 
    return graph;
	}

// Adding an edge to graph.
void addEdge(struct Graph* graph, int src, int dest, double weight, int *edge_no)
	{
		weight=floor(weight);
    // Add an edge from src to dest. The new node is added to the begining of adjList.
    graph->edge_array[*edge_no].next = graph->array[src].head;
	graph->edge_array[*edge_no].weight = weight;
	graph->edge_array[*edge_no].adjNode = dest;
    graph->array[src].head = graph->edge_array+(*edge_no);
    // Also, add an edge from dest to src as graph is undirected.
    *edge_no+=1;
    graph->edge_array[*edge_no].next = graph->array[dest].head;
	graph->edge_array[*edge_no].weight = weight;
	graph->edge_array[*edge_no].adjNode = src;
    graph->array[dest].head = graph->edge_array+(*edge_no);
	*edge_no+=1;
	}

void printGraph(struct Graph* graph)
	{
    int v;
    for (v = 0; v < graph->V; v++)
    	{
        struct AdjListNode* curr_node = graph->array[v].head;
        printf("%d:", v);
        while (curr_node)
        	{
            printf(" %d/%d", curr_node->adjNode, curr_node->weight);
            curr_node = curr_node->next;
        	}
        printf("\n");
    	}
	}

void get_numVertex(FILE *mapFile, int* num_V)
	{
	fscanf(mapFile, "%d", num_V);
	}

void get_numEdge(FILE *mapFile, int* num_E)
	{
	fscanf(mapFile, "%d", num_E);
	}

void read_coordinates(FILE *mapFile, int* num_V, short int* coor)
	{
	int vertex, i;
	for(i=0; i<*num_V; i++)
		{
		fscanf(mapFile, "%d", &vertex);
		fscanf(mapFile, "%hu", coor+(2*i));
		fscanf(mapFile, "%hu", coor+(2*i)+1);
		}
	}

void print_coordinates(int* num_V, short int* coor) //debug only
	{
	int i;
	for(i=0; i<*num_V; i++)
		{
		printf("%d\t%d\t%d\n",i,*(coor+i+i),*(coor+i+i+1));
		}
	}

double find_dist(int x1, int y1, int x2, int y2)
	{
	double dist;
	int x2_minus_x1 = x2-x1;
	int y2_minus_y1 = y2-y1;
	double squared_dist = (x2_minus_x1*x2_minus_x1)+(y2_minus_y1*y2_minus_y1);
	dist = sqrt(squared_dist);
	return dist;
	}

void build_adjList(FILE *mapFile, int* num_E, struct Graph* graph, short int* coor, int* edge_no)
	{
	int i;
	int src, dest;
	double weight=0.0;
	for(i=0; i<*num_E; i++)
		{
		fscanf(mapFile, "%d", &src);
		fscanf(mapFile, "%d", &dest);
		weight=find_dist(*(coor+(2*src)),*(coor+(2*src)+1),*(coor+(2*dest)),*(coor+(2*dest)+1));
    	addEdge(graph, src, dest, weight, edge_no);
		}
	}

// Heap Functions:

void minHeapify(struct MinHeap* minHeap, int root, struct Graph* graph)
	{
    int smallest = root;  			// Initialize smallest as root
    int left = (root << 1) + 1;  	// left = 2*root + 1
    int right = (root << 1) + 2; 	// right = 2*root + 2

    if (left < minHeap->size && minHeap->array[left].value < minHeap->array[smallest].value)
		{
		smallest = left;
		}
    if (right < minHeap->size && minHeap->array[right].value < minHeap->array[smallest].value)
		{
        smallest = right;
		}
    if (smallest != root) 			// check if swapping of root is needed.
    	{
		int vertex_a = minHeap->array[smallest].vertex;			// Begin Swap
		int vertex_b = minHeap->array[root].vertex;
		graph->array[vertex_a].node_heap = root;
		graph->array[vertex_b].node_heap = smallest;
		struct heapnode tmp = minHeap->array[smallest];
		minHeap->array[smallest] = minHeap->array[root];
		minHeap->array[root] = tmp;								// End Swap
        minHeapify(minHeap, smallest, graph);		
    	}
	}

struct MinHeap* BuildHeap(int init_val, int init_vertex, struct heapnode* array, struct Graph* graph) //used only once for initial vertex.
	{
    struct MinHeap* minHeap = (struct MinHeap*) malloc(sizeof(struct MinHeap));
    minHeap->size = 1;   	// initialize size of heap
	minHeap->array = array; 	// Assign address of first element of array
    minHeap->array[0].value = init_val; 	// PQ contains only the initial vertex.
    minHeap->array[0].vertex = init_vertex; 	// PQ contains only the initial vertex.
	graph->array[init_vertex].node_heap = 0;
    return minHeap;
	}

int extractMinPQ(struct MinHeap* minHeap, struct Graph* graph, int* c_ex_min)
	{
	int min=minHeap->array[0].vertex;
	minHeap->array[0]=minHeap->array[minHeap->size-1];
    minHeap->size--;  								// Reduce size of PQ.
	minHeapify(minHeap, 0, graph);					// Heapify root after swapping
	*(c_ex_min)+=1;
	return min;
	}

void decreaseKey(struct MinHeap* minHeap, int vertex, int key, struct Graph* graph, int* c_dec_key)
	{
	int heap_node = graph->array[vertex].node_heap;
	if (key > minHeap->array[heap_node].value)
		{
		return;
		}
	else
		{
		int parent=(heap_node-1)/2;
		minHeap->array[heap_node].value = key;
		while (heap_node>0 && minHeap->array[parent].value > minHeap->array[heap_node].value)
			{
			int vertex_a = minHeap->array[heap_node].vertex; // Begin Swap
			int vertex_b = minHeap->array[parent].vertex;
			graph->array[vertex_a].node_heap = parent;
			graph->array[vertex_b].node_heap = heap_node;
			struct heapnode tmp = minHeap->array[heap_node];
			minHeap->array[heap_node] = minHeap->array[parent];
			minHeap->array[parent] = tmp;					// End Swap
			heap_node=parent;
			parent=(heap_node-1)/2;
			}
		}*(c_dec_key)+=1;
	}

void insertPQ(struct MinHeap* minHeap, struct heapnode* array, int key, struct Graph* graph, int V, int* c_dec_key, int weight)
	{
	minHeap->size++;
	minHeap->array[minHeap->size-1].vertex = V;
	minHeap->array[minHeap->size-1].value = INT_MAX;
	graph->array[V].node_heap = minHeap->size-1;
	graph->array[V].distance = weight;
	decreaseKey(minHeap, V, key, graph, c_dec_key);
	}

void shortest_path(struct Graph* graph, int V, struct heapnode* Q, int src, int dest, short int* coor, int* c_ex_min,  int* c_dec_key)
	{
	int i=0;
	for(i=0; i<V; i++) 						// Initialize Single Source. 
		{
		graph->array[i].node_heap = -1;
		graph->array[i].visited = 0;
		}

	Q[src].value = 0; 						// make dist of source vertex 0.
	graph->array[src].distance = 0; 		// make dist of source vertex 0.
	graph->array[src].parent = src; 		// make parent of source = source.

	struct MinHeap* PQ = BuildHeap(0, src, Q, graph); 	// Build Priority Queue of all the vertices.
	int u = INT_MAX;
	while(PQ->size != 0 && u != dest ) 				
		{
		u = extractMinPQ(PQ, graph, c_ex_min);	// pick the min distance vertex from PQ.
		int vert = u;
		graph->array[u].visited=1;
    	struct AdjListNode* currNode = graph->array[u].head;
    	while (currNode)
			{vert = currNode->adjNode;
			double dist_vz = find_dist(*(coor+(2*vert)),*(coor+(2*vert)+1),*(coor+(2*dest)),*(coor+(2*dest)+1));
			double rem=fmod(dist_vz,1.0);
			int curr_wt; // we have to round-off curr_wt instead of using floor, because Astar algo gives a higher weight path compared to dijkstra's for few cases. This is because using the floor of heuristic distance, Astar might select the wrong path as its looking at a fake(lower) distance of incorrect path than what the actual path would have. Thus Astar returns a higher weight path corresponding to the incorrect path.
			if(rem>0.5)
			{curr_wt = graph->array[u].distance + currNode->weight + floor(dist_vz)+1;}
			else
			{curr_wt = graph->array[u].distance + currNode->weight + floor(dist_vz);}
			if(graph->array[vert].visited==0 && graph->array[vert].node_heap!=-1) // each enqueued adjacent vertex v.
				{
				if(PQ->array[graph->array[vert].node_heap].value > curr_wt)
					{
					decreaseKey(PQ, vert, curr_wt, graph, c_dec_key); // update the wt of currNode
					graph->array[vert].parent = u;
					graph->array[vert].distance=graph->array[u].distance + currNode->weight;
					}
				}
			if(graph->array[vert].visited==0 && graph->array[vert].node_heap==-1) // each remaining unenqueued adjacent vertex v
				{
				insertPQ(PQ, Q, curr_wt, graph, vert, c_dec_key, currNode->weight); // insert currNode in PQ with curr weight
				graph->array[vert].parent = u;
				graph->array[vert].distance=graph->array[u].distance + currNode->weight;
				}
				currNode = currNode->next;
			}
		}
	free(PQ);
	}

void printPath(struct Graph* graph, int dest, int src) // debug only
	{
int parent=graph->array[dest].parent;
	if(parent == src )
		{
		printf("%d ", src);
		return;
		}
	else
		{
		printPath(graph, parent, src);
		printf("%d ", parent);
		}
	}



int main(int argc, char *argv[])
{
	if (argc != 3) 
		{
		fprintf(stderr, "ERROR! Insufficient Inputs :( \n");
		fprintf(stderr, "Usage: ./SP_Astar.c <Map File> <query File>\n");
      	return EXIT_FAILURE;
      	}
	clock_t Start = 0;
	clock_t End = 0;
	double Read_Time = 0;
	Start = clock();

	FILE *mapFile, *queryFile;
	mapFile=fopen(argv[1],"r");
	if(mapFile == NULL) 
		{
		fprintf(stderr, "\n Error opening Map file %s \n", argv[1]);
		return EXIT_FAILURE;
		}
	queryFile=fopen(argv[2],"r");
	if(queryFile == NULL) 
		{
		fprintf(stderr, "\n Error opening query file %s \n", argv[2]);
		return EXIT_FAILURE;
		}
	int num_vertex; // # vertices
	int num_edge;   // # edges
	get_numVertex(mapFile, &num_vertex); // get number of vertices from map file.
	get_numEdge(mapFile, &num_edge); 	 // get number of edges from map file.
    struct Graph* graph = createGraph(num_vertex, num_edge); // create graph of num_vertex vertices.
	//printf("Number of vertices:%d, Edges:%d \n", num_vertex, num_edge);
	short int *coordinates = (short int*)malloc(2*num_vertex*sizeof(short int)); // to save coordinates. 
	read_coordinates(mapFile, &num_vertex, coordinates); 	// read coor from map file.
	//print_coordinates(&num_vertex, coordinates);
	int edge_ptr=0;											// ptr for edge_array in graph.
	build_adjList(mapFile, &num_edge, graph, coordinates, &edge_ptr); // create adj List.
    //printGraph(graph);

	struct heapnode *Q = (struct heapnode*)malloc(num_vertex*sizeof(struct heapnode)); // Priority Queue.
	int num_query;
	int src;
	int dest;
	fscanf(queryFile, "%d", &num_query);
	int i;
	int count_ex_min=0, count_dec_key=0;
	for(i=0; i<num_query; i++) // calculate distance for all queries.
		{
		fscanf(queryFile, "%d", &src);
		fscanf(queryFile, "%d", &dest);
		shortest_path(graph, num_vertex, Q, src, dest, coordinates, &count_ex_min, &count_dec_key);
		if(graph->array[dest].distance == INT_MAX || dest >= num_vertex || src >= num_vertex)
			{
			printf("INF \n");
			printf("%d %d \n",src,dest);
			}
		else
			{
			printf("%d\n",graph->array[dest].distance);
			printPath(graph,dest,src);
			printf("%d\n",dest);
			}
		}
	printf("Total extract min: %d \nTotal dec key: %d\n", count_ex_min, count_dec_key);

	fclose(mapFile);
	fclose(queryFile);
	free(coordinates);
	free(Q);
	free(graph->array);
	free(graph->edge_array);
	free(graph);
	End = clock();
	Read_Time += difftime(End, Start);
	printf("Query time: %le\n", Read_Time/CLOCKS_PER_SEC);
   	return EXIT_SUCCESS;
}

