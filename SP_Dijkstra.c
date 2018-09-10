#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <stdbool.h>

/******************************
*[Piyush Agarwal]
*[SP_Dijkstra.c]: [7/22/18].
*[agarw145@purdue.edu]
*Compiling: [gcc -Werror -lm -Wall -O3 SP_Dijkstra.c -o Dijkstra]
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
    struct AdjListNode *head; 
	};

struct Graph 			// structure for graph
	{
    int V;
    struct AdjList* array;
	};

// Create new adj list node.
struct AdjListNode* newAdjListNode(int adjNode)
	{
    struct AdjListNode* newNode = (struct AdjListNode*) malloc(sizeof(struct AdjListNode));
    newNode->adjNode = adjNode;
    newNode->next = NULL;
    return newNode;
	}

// Create graph of V vertices.
struct Graph* createGraph(int V)
	{
    struct Graph* graph = (struct Graph*) malloc(sizeof(struct Graph));
    graph->V = V;
    graph->array = (struct AdjList*) malloc(V * sizeof(struct AdjList));
 
    int i;
    for (i = 0; i < V; i++)
        {
		graph->array[i].parent = INT_MAX;
		graph->array[i].distance = INT_MAX;
		graph->array[i].node_heap = i;
		graph->array[i].head = NULL;
		}
 
    return graph;
	}

// Adding an edge to graph.
void addEdge(struct Graph* graph, int src, int dest, double weight)
	{double rem;
	rem=fmod(weight,1.0); // taking the floor of distance.
	//if(rem>0.5) // use it if you want to round off the distance to nearest integer.
	//	weight=floor(weight);
	//else
		weight=floor(weight);
    // Add an edge from src to dest. The new node is added to the begining of adjList.
    struct AdjListNode* newNode = newAdjListNode(dest);
    newNode->next = graph->array[src].head;
	newNode->weight = weight;
    graph->array[src].head = newNode;
    // Also, add an edge from dest to src as graph is undirected.
    newNode = newAdjListNode(src);
    newNode->next = graph->array[dest].head;
	newNode->weight = weight;
    graph->array[dest].head = newNode;
	}

void free_mem(struct Graph* graph, int num_V)
	{
	int v;
	for(v=0;v<num_V; v++)
		{
		free(graph->array[v].head);
		}
	}	

void printGraph(struct Graph* graph) // debug only
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
	int vertex, x_coor, y_coor, i;
	for(i=0; i<*num_V; i++)
		{
		fscanf(mapFile, "%d", &vertex);
		fscanf(mapFile, "%d", &x_coor);
		fscanf(mapFile, "%d", &y_coor);
		*(coor+(2*i))=x_coor;
		*(coor+(2*i)+1)=y_coor;
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

void build_adjList(FILE *mapFile, int* num_E, struct Graph* graph, short int* coor)
	{
	int i;
	int src, dest;
	double weight=0.0;
	for(i=0; i<*num_E; i++)
		{
		fscanf(mapFile, "%d", &src);
		fscanf(mapFile, "%d", &dest);
		weight=find_dist(*(coor+(2*src)),*(coor+(2*src)+1),*(coor+(2*dest)),*(coor+(2*dest)+1));
    	addEdge(graph, src, dest, weight);
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

struct MinHeap* BuildHeap(struct heapnode* array, int num_V, struct Graph* graph)
{
    int i;
    struct MinHeap* minHeap = (struct MinHeap*) malloc(sizeof(struct MinHeap));
    minHeap->size = num_V;   	// initialize size of heap
    minHeap->array = array; 	// Assign address of first element of array

    for (i = (minHeap->size - 2) / 2; i >= 0; i--)
        {
		minHeapify(minHeap, i, graph);
		}
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

void shortest_path(struct Graph* graph, int V, struct heapnode* Q, int src, int dest, int* c_ex_min, int* c_dec_key)
	{
	int i;
	for(i=0; i<V; i++) 						// Initialize Single Source. 
		{
		Q[i].value = INT_MAX; 				// v.d=INF
		Q[i].vertex = i;					// graph vertex to maintain handle to heap node.
		graph->array[i].parent = INT_MAX;
		graph->array[i].distance = INT_MAX;
		graph->array[i].node_heap = i;
		}

	Q[src].value = 0; 						// make dist of source vertex 0.
	graph->array[src].distance = 0; 		// make dist of source vertex 0.
	graph->array[src].parent = src; 		// make parent of source = source.


	struct MinHeap* PQ = BuildHeap(Q, V, graph); 	// Build Priority Queue of all the vertices.
	int u = INT_MAX;
	while(PQ->size != 0 && u != dest ) 				
		{
		u = extractMinPQ(PQ, graph, c_ex_min);	// pick the min distance vertex from PQ.
    	struct AdjListNode* currNode = graph->array[u].head;
    	while (currNode)
    		{
			int curr_vertex = currNode->adjNode;
			int curr_wt = currNode->weight;
			if (graph->array[u].distance+curr_wt < graph->array[curr_vertex].distance) // Relax vertices.
				{
				graph->array[curr_vertex].distance = graph->array[u].distance+curr_wt;
				graph->array[curr_vertex].parent = u;
				decreaseKey(PQ, curr_vertex, graph->array[u].distance+curr_wt, graph, c_dec_key);
				}
        	currNode = currNode->next;
    		}
		}
	}

void printPath(struct Graph* graph, int dest, int src)
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
		fprintf(stderr, "Usage: ./SP_Dijkstra.c <Map File> <query File>\n");
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
	get_numEdge(mapFile, &num_edge); // get number of edges from map file.
    struct Graph* graph = createGraph(num_vertex); // create graph of num_vertex vertices.
	short int *coordinates = (short int*)malloc(2*num_vertex*sizeof(short int)); // to save coordinates. 
	read_coordinates(mapFile, &num_vertex, coordinates); // read coor from map file.
	build_adjList(mapFile, &num_edge, graph, coordinates); // create adj List.
	free(coordinates);
    //printGraph(graph); //debug

	struct heapnode *Q = (struct heapnode*)malloc(num_vertex*sizeof(struct heapnode)); // Priority Queue.
	int num_query;
	int src;
	int dest;
	fscanf(queryFile, "%d", &num_query);
	int i;
	int count_ex_min=0, count_dec_key=0;
	for(i=0; i<num_query; i++)
		{
		fscanf(queryFile, "%d", &src);
		fscanf(queryFile, "%d", &dest);
		shortest_path(graph, num_vertex, Q, src, dest, &count_ex_min, &count_dec_key);
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
	free(Q);
	free_mem(graph, num_vertex);
	free(graph->array);
	free(graph);
	End = clock();
	Read_Time += difftime(End, Start);
	printf("Map Read time: %le\n", Read_Time/CLOCKS_PER_SEC);
   	return EXIT_SUCCESS;
}

