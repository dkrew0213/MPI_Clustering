/*
 * DNAStrand.h
 *
 *  Created on: Jul 29, 2014
 *      Author: Douglas Rew
 */

#ifndef DNASTRAND_H_
#define DNASTRAND_H_

#include <ctime>
#include <string>
#include <time.h>
using namespace std;

class DNAStrandHelper {
public:
	// This will be the constructor
	// Length will be the length of the DNA Strand
	// dataSize will be the size of the Data
	// clusters will be the number for clustering
	DNAStrandHelper();

	// Send the parms to every one.
	void setParams(int length, int dataSize, int centroidNum, double threshold ,int iter);

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

	// This will validate the datapoint with the centroids
	int validateStrandForCluster(char *strand, char **centriods);

	// This will recalculate the centroids
	void validateCentriod(char **data, char **centriods, int *cluster);

	// This will print out the data with the corresponding clusters
	void printAllWithCluster(char **data, int *cluster, string fileName);



	// MPI related functions
	// This will init the MPI
	void init_MPI();

	// This will start the classifying the data point into clusters
	void MPI_findCluster(char **data, char **centriods, int *cluster);

	// This will distribute the centroids
	void MPI_distributeCentroids();

	// This will send the data to all slaves
	void MPI_sendParitionData(char **data);

	// This will be used to receive data from RANK ==0
	void MPI_recievePartitionData();

	// This will close up the MPI
	void final_MPI();

	// This is used to partition the data size for the clusters
	int* partitionCluster();

	// This will return the rank
	int MPI_getRank();

	// This is used to send the classification information
	void MPI_sendClusters();

	// This is used to recieve the classification information
	void MPI_recieveClusters(int *cluster);

private:
	// Length for the DNA strand
	int length;
	// Total size of the data
	int dataSize;
	// Number of the centroid
	int centroidNum;

	// for MPI
	// My rank
	int myRank;
	// The size of the whole world in MPI
	int myWorld;
	// This will be my local data size
	int partitionedDataSize;
	// This will be the partition data
	char **partitionedData;
	// This will be the local centroids
	char **centriods;
	// This will be the local classification information
	int *partitionedCluster;

	// Local copy for the threshold
	double threshold;
	// Local copy for the max iteration
	int iteration;

};

#endif /* DNASTRAND_H_ */
