/*
* Code by Aaroh Sharma and Kamil Yildirim
* 
* EIDs: as225925, ky5637
*/


#include <iostream>
#include <fstream>
#include <vector>
#include <float.h>
#include <algorithm>
#include <boost/format.hpp>
#include <boost/histogram.hpp>

#define DAMPING_FACTOR 0.85
#define EPSILON 0.0001
#define WIKI_PATH "data/wiki.dimacs"
#define ROAD_NY_PATH "data/road-NY.dimacs"
#define RMAT15_PATH "data/rmat15.dimacs"

using namespace std;
// using namespace boost::histogram;

struct Graph {
    vector<double> labels;                   // List of node labels
    vector<int> offsets;                  // List of to pointers
    vector<pair<int, int>> dest_weights;  // List of dest and weight pairs
};

Graph COOtoCSR (vector<vector<int>> cooRep, int numNodes, int numEdges);
Graph readGraphFromDIMACS(const string& filePath);
double pushPageRankIteration(Graph& graph, double dampingFactor);

bool sortByFrom(const vector<int>& a, const vector<int>& b) {
    return a[0] < b[0];
}

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

    vector<vector<int>> cooRep(numEdges, vector<int>(3, 0));
    
    // Add edges to adjacency matrix
    int count = 0;
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
            cooRep[count++] = {from, to, weight};
        }
    }

    sort(cooRep.begin(), cooRep.end(), sortByFrom);
    Graph graph = COOtoCSR(cooRep, numNodes, numEdges);

    file.close();
    return graph;
}

// Convert COO representation to CSR representation
Graph COOtoCSR (vector<vector<int>> cooRep, int numNodes, int numEdges) {
    Graph graph;
    graph.offsets.push_back(0);
    graph.labels.push_back(0.0);

    double initPageRank = 1 / (double) numNodes;

    // Iterate over the COO representation and add the edges to the CSR representation
    int curIndex = 0;
    for (int i = 1; i <= numNodes; i++) {
        graph.labels.push_back(initPageRank);
        graph.offsets.push_back(graph.dest_weights.size());
        for (; curIndex < numEdges; curIndex++) {
            if (cooRep[curIndex][0] != i) {
                break;
            }
            graph.dest_weights.push_back(make_pair(cooRep[curIndex][1], cooRep[curIndex][2]));
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

    file << "p sp " << numVertices << " " << numEdges << endl;

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
    double sum = 0.0;

    for (int i = 1; i < numVertices + 1; i++) {
        sum += g.labels[i];
        file << i << " " << g.labels[i] << endl; // Node numbers start from 1
    }

    file.close();
}

void pageRank(Graph& graph, double dampingFactor, double epsilon) {
    double prevMaxDiff = DBL_MAX;
    while (prevMaxDiff > epsilon) {
        prevMaxDiff = pushPageRankIteration(graph, dampingFactor);
    }

    double totalRank = 0;
    int numVertices = graph.offsets.size() - 2;

    for (int i = 1; i < numVertices + 1; i++) {
        totalRank += graph.labels[i];
    }
    // Normalize the PageRank scores
    for (int i = 1; i < numVertices + 1; i++) {
        graph.labels[i] = graph.labels[i] / totalRank;
    }
}

// Function to perform one iteration of the push-style PageRank algorithm
double pushPageRankIteration(Graph& graph, double dampingFactor) {
    int numVertices = graph.offsets.size() - 2;
    vector<double> newPageRank(numVertices + 1, 0.0);

    for (int i = 1; i < numVertices + 1; i++) {
        int outDegree = graph.offsets[i + 1] - graph.offsets[i];

        if (outDegree > 0) {
            double pushWeight = graph.labels[i];
            double totalWeight = 0.0;
            for (int j = graph.offsets[i]; j < graph.offsets[i + 1]; j++) {
                totalWeight += graph.dest_weights[j].second;
            }
            
            if (totalWeight > 0)
                pushWeight = pushWeight / totalWeight;

            for (int j = graph.offsets[i]; j < graph.offsets[i + 1]; j++) {
                int dest = graph.dest_weights[j].first;
                int weight = graph.dest_weights[j].second;
                if (weight > 0)
                    newPageRank[dest] += pushWeight * weight;
            }
        }
    }

    // Apply damping factor and update PageRank
    double teleportWeight = (1.0 - dampingFactor) / numVertices;
    double maxDiff = -DBL_MAX;

    for (int i = 1; i < numVertices + 1; i++) {
        newPageRank[i] = dampingFactor * newPageRank[i] + teleportWeight;
    }

    // Normalize the PageRank scores
    for (int i = 1; i < numVertices + 1; i++) {
        maxDiff = max(maxDiff, abs(newPageRank[i] - graph.labels[i]));
    }

    graph.labels = newPageRank;
    return maxDiff;
}

// void writeHistogramToFile (Graph& g, const string& filename) {
//     ofstream file(filename);
//     if (!file.is_open()) {
//         cerr << "Error: Unable to open the file: " << filename << endl;
//         return;
//     }

//     int numVertices = g.offsets.size() - 2;
//     auto h = make_histogram(axis::integer<>(1, numVertices + 1, "Out edges connected to Nodes"));

//     for (int i = 1; i < numVertices + 1; i++) {
//         int outDegree = g.offsets[i + 1] - g.offsets[i];
//         h(outDegree);
//     }

//     for (auto&& x : indexed(h)) {
//         file << boost::format("out-degree of %i: %i\n") % (x.index() + 1) % *x;
//     }
// }


// Function to compute the histogram of outgoing edges for the graph
void computeOutgoingEdgesHistogram(const Graph& graph, const std::string& filename) {
    ofstream file(filename);
    if (!file.is_open()) {
        cerr << "Error: Unable to open the file: " << filename << endl;
        return;
    }

    int numVertices = graph.offsets.size() - 2;
    std::vector<int> out_degrees(numVertices + 1, 0);

    for (int i = 1; i < graph.offsets.size() - 1; i++) {
        int outgoingEdges = graph.offsets[i + 1] - graph.offsets[i];
        out_degrees[i] = outgoingEdges;
    }

    int max_value = *max_element(out_degrees.begin(), out_degrees.end());
    int min_value = *min_element(out_degrees.begin(), out_degrees.end());

    int numBins = max_value - min_value + 1 < 10 ? max_value - min_value + 1 : 10;

    auto axis = boost::histogram::axis::regular<>(numBins, min_value, max_value + 1, "Outgoing Edges");
    auto h = boost::histogram::make_histogram(axis);

    // Fill the histogram with PageRank values
    for (int i = 1; i < out_degrees.size(); i++) {
        h(out_degrees[i]);
    }

    if (numBins == 10) {
        for (auto&& x : indexed(h)) {
            file << boost::format("Bin %i [%i, %i]: Count = %i\n") % x.index() 
                % (int) x.bin().lower() % (int) x.bin().upper() % *x;
        }
    } else {
        for (auto&& x : indexed(h)) {
            file << boost::format("Bin %i [%i]: Count = %i\n") % x.index() 
                % x.bin().lower() % *x;
        }
    }
}


// Function to create and show a histogram of PageRank values
// Function to create a PageRank histogram
void createPageRankHistogram(Graph& graph, int numBins, const std::string& filename) {
    ofstream file(filename);
    if (!file.is_open()) {
        cerr << "Error: Unable to open the file: " << filename << endl;
        return;
    }

    vector<double> pageRanks = graph.labels;
    double max_value = (*max_element(pageRanks.begin(), pageRanks.end())) * 1.001;
    double min_value = (*min_element(pageRanks.begin(), pageRanks.end())) * 0.009;
    auto axis = boost::histogram::axis::regular<>(numBins, min_value, max_value, "Page Ranks");
    auto h = boost::histogram::make_histogram(axis);

    // Fill the histogram with PageRank values
    for (int i = 1; i < pageRanks.size(); i++) {
        h(pageRanks[i]);
    }

    for (auto&& x : indexed(h)) {
        file << boost::format("Bin %i [%i, %i]: Count = %i\n") % x.index() 
            % x.bin().lower() % x.bin().upper() % *x;
    }
}

void runWiki() {
    Graph graph = readGraphFromDIMACS(WIKI_PATH);
    CSRtoDIMACS(graph, "outputs/wiki/wiki_copy.dimacs");
    printNodeNumbersAndLabels(graph, "outputs/wiki/init_labels.txt");

    pageRank(graph, DAMPING_FACTOR, EPSILON);
    printNodeNumbersAndLabels(graph, "outputs/wiki/final_labels.txt");

    createPageRankHistogram(graph, 10, "outputs/wiki/page_ranks_histogram.txt");
    computeOutgoingEdgesHistogram(graph, "outputs/wiki/outgoing_edges_histogram.txt");
}

void runRoadNY () {
    Graph graph = readGraphFromDIMACS(ROAD_NY_PATH);
    CSRtoDIMACS(graph, "outputs/road_NY/road_ny_copy.dimacs");
    printNodeNumbersAndLabels(graph, "outputs/road_NY/init_labels.txt");

    pageRank(graph, DAMPING_FACTOR, EPSILON);
    printNodeNumbersAndLabels(graph, "outputs/road_NY/final_labels.txt");

    createPageRankHistogram(graph, 10, "outputs/road_NY/page_ranks_histogram.txt");
    computeOutgoingEdgesHistogram(graph, "outputs/road_NY/outgoing_edges_histogram.txt");
}

void runRMAT15 () {
    Graph graph = readGraphFromDIMACS(RMAT15_PATH);
    CSRtoDIMACS(graph, "outputs/rmat15/rmat15_copy.dimacs");
    printNodeNumbersAndLabels(graph, "outputs/rmat15/init_labels.txt");

    pageRank(graph, DAMPING_FACTOR, EPSILON);
    printNodeNumbersAndLabels(graph, "outputs/rmat15/final_labels.txt");

    createPageRankHistogram(graph, 10, "outputs/rmat15/page_ranks_histogram.txt");
    computeOutgoingEdgesHistogram(graph, "outputs/rmat15/outgoing_edges_histogram.txt");
}

int main(int argc, char* argv[]) {
    runWiki();
    runRoadNY();
    runRMAT15();
    return 0;
}