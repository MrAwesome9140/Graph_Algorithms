To run the program:
1.First run "make all" in the command line.
2.Then run ./graph_algorithm to run the executable and get the output.
3.Then, "make clean" can be run to clean up and remove the binary.

Folder Description:

Each folder for each of the tests have a final_labels, init_labels, a copy.dimacs file, outgoing_edges_histogram, and page_ranks_histogram. 
The init_labels file hold the initial page rank values before the page rank algorithm runs. The final_labels holds the values after 
the page rank algorithm runs, hence the final labels. The copy.dimacs file for each of the tests hold the dimacs graph that 
we get after we convert the original dimacs to csr and then back to dimacs to print it out. The outgoing_edges_histogram holds the 
histogram of the number of outgoing edges connected to each node. The page_ranks_histogram holds the histogram of all of the page ranks
of each node. 
