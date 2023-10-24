#include <iostream>
#include <fstream>
#include <vector>
#include <float.h>
#include <algorithm>
#include <boost/format.hpp>
#include <boost/histogram.hpp>

#define DAMPING_FACTOR 0.85
#define EPSILON 0.0001

using namespace std;
using namespace boost::histogram;

struct Graph {
    vector<double> labels;                   // List of node labels
    vector<int> offsets;                  // List of to pointers
    vector<pair<int, int>> dest_weights;  // List of dest and weight pairs
};

Graph adjacencyToCSR (vector<vector<int>> adjacency, int numNodes, int numEdges);
Graph readGraphFromDIMACS(const string& filePath);
double pushPageRankIteration(Graph& graph, double dampingFactor);

// Function to read a graph in DIMACS format and construct CSR representation
Graph readGraphFromDIMACS(const string& filePath) {
    ifstream file(filePath);
    if (!file.is_open()) {
        cerr << "Error: Unable to open the file: " << filePath << endl;
        exit(-1);
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
            break;
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
    graph.offsets.push_back(0);
    graph.labels.push_back(0.0);

    double initPageRank = 1 / (double) numNodes;

    // Iterate over the adjacency matrix and add the edges to the CSR representation
    for (int i = 1; i < numNodes + 1; i++) {
        graph.labels.push_back(initPageRank);
        graph.offsets.push_back(graph.dest_weights.size());
        for (int j = 1; j < numNodes + 1; j++) {
            if (adjacency[i][j] != 0) {
                graph.dest_weights.push_back(make_pair(j, adjacency[i][j]));
            }
        }
    }
    graph.offsets.push_back(graph.dest_weights.size());

    return graph;
}

void CSRtoDIMACS (Graph g, const string& filename) {
    ofstream file(filename);
    if (!file.is_open()) {
        cerr << "Error: Unable to open the file: " << filename << endl;
        return;
    }

    int numVertices = g.offsets.size() - 2;
    int numEdges = g.dest_weights.size();

    file << "p sp " << numVertices << " " << numEdges << std::endl;

    for (int i = 1; i < numVertices + 1; i++) {
        for (int j = g.offsets[i]; j < g.offsets[i + 1]; j++) {
            file << "a " << i << " " << g.dest_weights[j].first << " " << g.dest_weights[j].second << std::endl;
        }
    }
}


// Function to print node numbers and labels to a file
void printNodeNumbersAndLabels(const Graph& g, const string& filename) {
    ofstream file(filename);
    if (!file.is_open()) {
        cerr << "Error: Unable to open the file: " << filename << endl;
        return;
    }

    int numVertices = g.offsets.size() - 2;

    for (int i = 1; i < numVertices + 1; i++) {
        file << i << " " << g.labels[i] << endl; // Node numbers start from 1
    }

    file.close();
}

void pageRank(Graph& graph, double dampingFactor, double epsilon) {
    double prevMaxDiff = DBL_MAX;
    while (prevMaxDiff > epsilon) {
        prevMaxDiff = pushPageRankIteration(graph, dampingFactor);
    }
}

// Function to perform one iteration of the push-style PageRank algorithm
double pushPageRankIteration(Graph& graph, double dampingFactor) {
    int numVertices = graph.offsets.size() - 2;
    vector<double> newPageRank(numVertices + 1, 0.0);

    for (int i = 1; i < numVertices + 1; i++) {
        int outDegree = graph.offsets[i + 1] - graph.offsets[i];

        if (outDegree > 0) {
            // double totalOutWeight = 0.0;
            // for (int j = graph.offsets[i]; j < graph.offsets[i + 1]; j++) {
            //     totalOutWeight += graph.dest_weights[j].second;
            // }

            // double pushWeight = graph.labels[i] / totalOutWeight;

            // for (int j = graph.offsets[i]; j < graph.offsets[i + 1]; j++) {
            //     int dest = graph.dest_weights[j].first;
            //     int weight = graph.dest_weights[j].second;
            //     newPageRank[dest] += pushWeight * weight;
            // }

            double pushWeight = graph.labels[i] / outDegree;

            for (int j = graph.offsets[i]; j < graph.offsets[i + 1]; j++) {
                int dest = graph.dest_weights[j].first;
                int weight = graph.dest_weights[j].second;
                newPageRank[dest] += pushWeight;
            }
        }
    }

    // Apply damping factor and update PageRank
    double teleportWeight = (1.0 - dampingFactor) / numVertices;
    double totalRank = 0;
    double maxDiff = -DBL_MAX;

    for (int i = 1; i < numVertices + 1; i++) {
        newPageRank[i] = dampingFactor * newPageRank[i] + teleportWeight;
        totalRank += newPageRank[i];
    }

    // Normalize the PageRank scores
    for (int i = 1; i < numVertices + 1; i++) {
        newPageRank[i] = newPageRank[i] / totalRank;
        maxDiff = max(maxDiff, abs(newPageRank[i] - graph.labels[i]));
    }

    graph.labels = newPageRank;
    return maxDiff;
}

void writeHistogramToFile (Graph& g, const string& filename) {
    ofstream file(filename);
    if (!file.is_open()) {
        cerr << "Error: Unable to open the file: " << filename << endl;
        return;
    }

    int numVertices = g.offsets.size() - 2;
    auto h = make_histogram(axis::integer<>(1, numVertices + 1, "Out edges connected to Nodes"));

    for (int i = 1; i < numVertices + 1; i++) {
        int outDegree = g.offsets[i + 1] - g.offsets[i];
        h(outDegree);
    }

    for (auto&& x : indexed(h)) {
        file << boost::format("out-degree of %i: %i\n") % (x.index() + 1) % *x;
    }
}


int main() {
    std::string filePath = "wiki.dimacs";
    Graph graph = readGraphFromDIMACS(filePath);
    CSRtoDIMACS(graph, "wiki2.dimacs");
    printNodeNumbersAndLabels(graph, "init_wiki_labels.txt");
    pageRank(graph, DAMPING_FACTOR, EPSILON);
    printNodeNumbersAndLabels(graph, "final_wiki_labels.txt");
    writeHistogramToFile(graph, "histogram.txt");
    return 0;
}