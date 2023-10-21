#include <iostream>
#include <fstream>
#include <vector>

using namespace std;

struct Graph {
    vector<int> pointers;               // List of to pointers
    vector<pair<int, int>> dest_weights;    // List of dest and weight pairs
};

// Function to read a graph in DIMACS format and construct CSR representation
Graph readGraphFromDIMACS(const string& filePath) {
    ifstream file(filePath);
    if (!file.is_open()) {
        cerr << "Error: Unable to open the file: " << filePath << endl;
        exit(1);
    }
    int numNodes, numEdges;

    // Get problem information, number of nodes and edges
    string line;
    while (getline(file, line)) {
        if (line.empty() || line[0] == 'c') {
            // Skip comments
            continue;
        }

        if (line[0] == 'p') {
            // Parse the problem line to get the number of nodes and edges
            sscanf(line.c_str(), "p sp %d %d", &numNodes, &numEdges);
        }
    }

    vector<vector<int>> adjacencyMatrix(numNodes + 1, 
        vector<int>(numNodes + 1, 0));
    
    // Add edges to adjacency matrix
    while (getline(file, line)) {
        if (line.empty() || line[0] == 'c') {
            // Skip comments
            continue;
        }

        if (line[0] == 'a') {
            // Parse an edge line
            int from, to, weight;
            std::sscanf(line.c_str(), "a %d %d %d", &from, &to, &weight);

            // Store the destination node in the adjacency list
            adjacencyMatrix[from][to] = weight;
        }
    }

    Graph graph = adjacencyToCSR(adjacencyMatrix, numNodes, numEdges);

    file.close();
    return graph;
}

/* Function to convert adjacency matrix to CSR representation
*/
Graph adjacencyToCSR (vector<vector<int>> adjacency, int numNodes, int numEdges) {
    Graph graph;
    graph.pointers.push_back(0);

    // Iterate over the adjacency matrix and add the edges to the CSR representation
    for (int i = 1; i < numNodes; i++) {
        for (int j = 1; j < numNodes; j++) {
            if (adjacency[i][j] != 0) {
                graph.dest_weights.push_back(
                    make_pair(j, adjacency[i][j])
                );
            }
        }
        graph.pointers.push_back(graph.dest_weights.size());
    }

    return graph;
}

int main() {
    std::string filePath = "wiki.dimacs";
    Graph graph = readGraphFromDIMACS(filePath);

    return 0;
}