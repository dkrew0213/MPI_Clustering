/*
 * DNAStrand.h
 *
 *  Created on: Jul 29, 2014
 *      Author: Douglas Rew
 */

#ifndef DNASTRAND_H_
#define DNASTRAND_H_

#include <string>
using namespace std;

class DNAStrandHelper {
public:
	// This will be the constructor
	// Length will be the length of the DNA Strand
	// dataSize will be the size of the Data
	// clusters will be the number for clustering
	DNAStrandHelper(int length, int dataSize, int clusters);

	// Destructor
	virtual ~DNAStrandHelper();


	// This will return a randomType for the DNA
	char getRandomType();

	// This will return a random DNA Strand
	char* getRamdomDNAStrand();

	// This will print out the DNA strand
	void printStrand(char* strand);

	// This will compare two strand
	int compareTwo(char* strandA, char* strandB);

	// This will pick random centriods
	char** pickCentriods(char** data);

	// This will find the clusters using the threshold and the maxiteration
	void findCluster(char **data, char **centriods, int *cluster, double threshold, int maxIteration);

	// This will validate the datapoint with the centroids
	int validateStrandForCluster(char *strand, char **centriods);

	// This will recalculate the centroids
	void validateCentriod(char **data, char **centriods, int *cluster);

	// This will print out the data with the corresponding clusters
	void printAllWithCluster(char **data, int *cluster, string fileName);

private:
	// Length for DNA Strand
	int length;
	// Total Data size
	int dataSize;
	// Number of clusters (centroids)
	int clusters;
};

#endif /* DNASTRAND_H_ */
