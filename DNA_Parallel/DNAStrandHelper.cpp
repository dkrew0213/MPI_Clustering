/*
 * DNAStrand.cpp
 *
 *  Created on: Jul 29, 2014
 *      Author: Douglas Rew
 */

#include "DNAStrandHelper.h"
#include <iostream>
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */
#include <math.h>
#include <limits>
#include "mpi.h"
#include <fstream>
#include <string>

using namespace std;

// Constructor
DNAStrandHelper::DNAStrandHelper() {

	srand(time(0));

	// for MPI
	this->myRank = 0;
	this->myWorld =0;
	this->partitionedDataSize = 0;
	this->partitionedData = NULL;
	this->centriods= NULL;
	this->partitionedCluster = NULL;

	this->length = 0;
	this->dataSize = 0;
	this->centroidNum = 0;



}

// Function that will send the parameter to all the child from the master which is the RANK ==0
void DNAStrandHelper::setParams(int length, int dataSize, int centroidNum,  double threshold , int iter){

	if(this->myRank == 0){
		this->length = length;
		this->dataSize = dataSize;
		this->centroidNum = centroidNum;
		this->threshold = threshold;
		this->iteration = iter;

		for(int i = 1; i< this->myWorld; i++){
			 MPI_Send(&this->centroidNum, 1, MPI_INT, i, 0,MPI_COMM_WORLD);
			 MPI_Send(&this->dataSize, 1, MPI_INT, i, 0,MPI_COMM_WORLD);
			 MPI_Send(&this->length, 1, MPI_INT, i, 0,MPI_COMM_WORLD);
			 MPI_Send(&this->threshold, 1, MPI_DOUBLE, i, 0,MPI_COMM_WORLD);
			 MPI_Send(&this->iteration, 1, MPI_INT, i, 0,MPI_COMM_WORLD);
		}

	} else {
		  MPI_Recv(&this->centroidNum, 1, MPI_INT, 0, 0, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
		  MPI_Recv(&this->dataSize, 1, MPI_INT, 0, 0, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
		  MPI_Recv(&this->length, 1, MPI_INT, 0, 0, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
		  MPI_Recv(&this->threshold, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
		  MPI_Recv(&this->iteration, 1, MPI_INT, 0, 0, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
	}

}

// Destructor
// For the parallelization we are carrying the partition data into its own instant
// we are cleaning this up here.
DNAStrandHelper::~DNAStrandHelper() {


	if( this->partitionedData != NULL){
		for(int i = 0; i < this->partitionedDataSize ; i++){
			delete[] partitionedData[i];
		}
		delete[] partitionedData;
	}
	if(this -> centriods != NULL){
		for(int i = 0; i < this->centroidNum; i++){
			delete[] centriods[i];
		}
	}
	if( this->partitionedCluster != NULL){
		delete[] partitionedCluster;
	}
}

// This will generate the DNA Strand
char* DNAStrandHelper::getRamdomDNAStrand(){
	char* strand = new char[this->length];
	for(int i = 0 ; i < this->length; i++){
		strand[i] = this->getRandomType();
	}
	return strand;
}

// This will genreate a random DNA alphabet
char DNAStrandHelper::getRandomType(){

	int index = rand()%4;
	switch(index){
	case 0:
		return 'A';
	case 1:
		return 'G';
	case 2:
		return 'C';
	case 3:
		return 'T';
	}
	return '-';

}

// Helper function to print out the strand
void DNAStrandHelper::printStrand(char* strand){
	for(int i = 0; i < this->length; i++){
		char value = strand[i];
		switch(value){
		case 'A':
			cout<< "A"; break;
		case 'G':
			cout <<"G"; break;
		case 'C':
			cout<<"C"; break;
		case 'T':
			cout<<"T"; break;
		default:
			cout<< '*'; break;
		}
	}

	cout<<endl;
}

// Function that will compare the two DNA Strand
int DNAStrandHelper::compareTwo(char* strandA, char* strandB){
	int counter = 0;

	for(int i = 0; i < this->length; i++){
		if(strandA[i] != strandB[i]){
			counter++;
		}
	}

	return counter;
}

// This will pick random centroids from the data points
char** DNAStrandHelper::pickCentriods(char** data){

	char** centroids = new char*[this->centroidNum]();

	for(int i = 0; i < this->centroidNum; i++){
		centroids[i] = new char[this->length]();

		int index = rand()%this->dataSize;
		for(int j = 0; j < this->length; j++){
			centroids[i][j] = data[index][j];
		}
	}
	return centroids;
}




// This will validate the DNA strand and determine which cluster it belongs to.
int DNAStrandHelper::validateStrandForCluster(char *strand, char **centriods){
	int index = 0;
	int minDistance = this->length * 2;
	for(int i = 0; i < this->centroidNum; i++){
		// getting the distance of the two
		int distance = this->compareTwo(centriods[i], strand);
		if( distance < minDistance){
			index = i;
			minDistance = distance;
		}

	}
	return index;
}

// This will recalculate the centroids
void DNAStrandHelper::validateCentriod(char **data, char **centriods, int *cluster){

	// Setting up a middle object so that we could so the counting easily
	double ***middleStrand = new double**[this->centroidNum];
	for(int mid = 0; mid < this->centroidNum; mid++){
		middleStrand[mid] = new double*[this->length];
		for(int len = 0; len<this->length; len++){
			for(int c = 0; c < 4; c++){
				middleStrand[mid][len] = new double[4];
				middleStrand[mid][len][0] = 0;
				middleStrand[mid][len][1] = 0;
				middleStrand[mid][len][2] = 0;
				middleStrand[mid][len][3] = 0;

			}

		}
	}

	// Keep up with the cluster count
	int *clusterCounter = new int[this->centroidNum];
	for(int index = 0; index < this->centroidNum; index++){
		clusterCounter[index] = 0;
	}

	for(int i = 0; i < this->centroidNum; i++){
		// computing the per centroid
		for(int j = 0; j < this->dataSize; j++){
			if(cluster[j] == i){
				clusterCounter[i] = clusterCounter[i]+1;
				for(int k = 0; k < this->length; k++){
					if( data[j][k] == 'A'){ // taking it highest count
						middleStrand[i][k][0] = middleStrand[i][k][0] + 1;
					} else if( data[j][k] == 'G'){
						middleStrand[i][k][1] = middleStrand[i][k][1] + 1;
					} else if( data[j][k] == 'C'){
						middleStrand[i][k][2] = middleStrand[i][k][2] + 1;
					} else if( data[j][k] == 'T'){
						middleStrand[i][k][3] = middleStrand[i][k][3] + 1;
					}
				}
			}
		}
	}

	// Picking the most common value in the middleStrand
	for(int i = 0; i < this->centroidNum; i++){
		for(int k = 0; k < this->length; k++){

			int index = 0;
			int minCount = numeric_limits<int>::min();
			for(int c = 0; c < 4 ; c++){
				if(minCount < middleStrand[i][k][c] ){
					minCount = middleStrand[i][k][c];
					index = c;
				}
			}


			if(index == 0){
				centriods[i][k] = 'A';
			} else if(index == 1){
				centriods[i][k] = 'G';
			} else if(index == 2){
				centriods[i][k] = 'C';
			} else if(index == 3){
				centriods[i][k] = 'T';
			} else {
				centriods[i][k] = '*';
			}

		}

	}


	// Cleanin up the middle data
	for(int mid = 0; mid < this->centroidNum; mid++){
		for(int len = 0; len<this->length; len++){
				delete[] middleStrand[mid][len];
			}
			delete[] middleStrand[mid];
	}


	delete[] middleStrand;
	delete[] clusterCounter;
}

// Saving all data into a file.
void DNAStrandHelper::printAllWithCluster(char **data, int *cluster, string fileName){

	  ofstream myfile;
	  myfile.open(fileName.c_str(), ios::out | ios::app | ios::binary );
	  if(myfile.is_open()){
			for(int i =0 ; i < this->dataSize; i++){
				myfile << "Cluster : " << cluster[i] << "  = ";
				for(int j = 0; j < this->length; j++){
					char value = data[i][j];
					switch(value){
					case 'A':
						myfile<< "A"; break;
					case 'G':
						myfile <<"G"; break;
					case 'C':
						myfile<<"C"; break;
					case 'T':
						myfile<<"T"; break;
					default:
						myfile<< value; break;
					}
				}

				myfile<<endl;
			}
			myfile.close();
	  } else {
		  cout << "Error in opening the file = " << fileName << endl;
	  }

}


// This function wil distribute the centroids from the RANK == 0 to all others.
void DNAStrandHelper::MPI_distributeCentroids(){
	 if(this->myRank == 0){
		 for(int i = 1 ; i < this->myWorld; i++){
			 for(int j = 0; j < this->centroidNum; j++){
				 MPI_Send(this->centriods[j], this->length, MPI_CHAR, i, 0,MPI_COMM_WORLD);
			 }
		 }
	  } else {
		  if(this->centriods == NULL){
			  this->centriods = new char*[this->centroidNum];
			  for(int i = 0; i < this->centroidNum; i++){
				  this->centriods[i] = new char[this->length];
			  }
		  }

		  for(int j = 0; j < this->centroidNum; j++){;
			  MPI_Recv(this->centriods[j], this->length, MPI_CHAR, 0, 0, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
		  }
	  }

}


// This is where we classify the data point into a cluster
void DNAStrandHelper::MPI_findCluster(char **data, char **centriods, int *cluster){



	int iteration = 0;

	// RANK ==0 picks centroids from the whole data set.
	// How ever we want to use the same structure through out all instants
	// So we are copying a centroids into the local centroids.
	if(this->myRank == 0){

		this->centriods = new char*[this->centroidNum];
		for(int i =0; i < this->centroidNum; i++){
			this->centriods[i] = new char[this->length];
		}

		for(int i =0; i < this->centroidNum; i++){
			for(int j =0; j < this->length; j++){
				this->centriods[i][j] = centriods[i][j];
			}

		}
	}

	// Send out the centroids to people
	this->MPI_distributeCentroids();

	// This is a one time sending of the data
	if(this->myRank == 0){
		this->MPI_sendParitionData(data);
		cout << "Data all Sent" << endl;
	} else {
		this->MPI_recievePartitionData();
		cout << "Process:" << this->myRank << " Data all received " << endl;
	}

	// Initializing the partition cluster. This will contain an array of the classificatin of the clusters.
	if(this->partitionedDataSize != 0 ){
		this->partitionedCluster = new int[this->partitionedDataSize];
	}

	double clusterUpdated = this->partitionedDataSize;
	double comparedValue = this->dataSize * this->threshold; // We are using the total data size for doing threshold
	while(   clusterUpdated > comparedValue  && this->iteration > iteration ){
		clusterUpdated = 0;

		for(int i = 0; i < this->partitionedDataSize; i++){
				// this will give a validated cluster for the data point
			int validatedCluster = this->validateStrandForCluster(this->partitionedData[i], this->centriods);
			if( this->partitionedCluster[i] != validatedCluster){
				this->partitionedCluster[i] = validatedCluster;
				clusterUpdated++;
			}

		}

		iteration++;



		// Receive and send out the classification information
		if(this->myRank == 0){
			this->MPI_recieveClusters(cluster);
		} else {
			this->MPI_sendClusters();
		}

		if(this->myRank == 0){
			// The re calculation of the centroid happens in the RANK == 0;
			this->validateCentriod(data,this->centriods,cluster);
		}
		// We distribute the new centroids
		this->MPI_distributeCentroids();


		//Updating all with the clusterUpdated value
		// By doing this we will be able to keep up with handling the threshold
		MPI_Allreduce(&clusterUpdated,&clusterUpdated, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

		if(this->myRank == 0){
			// informational print
			cout <<"Info from RANK = 0 ::::::: Current Changing Count = " << clusterUpdated << " Threshold = " <<  comparedValue << endl;
		}
	}
}


// Setting up MPI
void DNAStrandHelper::init_MPI(){

	MPI_Init(NULL, NULL);
	MPI_Comm_rank(MPI_COMM_WORLD, &this->myRank);
	MPI_Comm_size(MPI_COMM_WORLD, &this->myWorld);
}

// Closing the MPI
void  DNAStrandHelper::final_MPI(){
	 MPI_Finalize();
}


// Funtion that will send partition data to all slave
void DNAStrandHelper::MPI_sendParitionData(char **data){
	cout << "Sending data to the nodes."<< endl;

	int *partitioned = this->partitionCluster();

	int dataIndex = 0;
	for(int i = 0; i < this->myWorld; i++){
		if(i == 0){

			// Initializing the internal partition data
			this->partitionedDataSize = partitioned[i];
			this->partitionedData = new char*[this->partitionedDataSize];

			for(int j = 0; j < this->partitionedDataSize; j++){
				this->partitionedData[j] = new char[this->length];
				for(int k = 0; k < this->length; k++){
					this->partitionedData[j][k] = data[dataIndex][k];
				}
				dataIndex++;
			}
		} else {

			int parition = partitioned[i];
			// We send the partition size first
			MPI_Send(&parition, 1, MPI_INT, i, 0,MPI_COMM_WORLD);

			// Sending the actual data
			for(int j = 0; j < parition; j++  ){
				MPI_Send(data[dataIndex], this->length, MPI_CHAR, i, 0,MPI_COMM_WORLD);
				dataIndex++;
			}
		}


	}

}


// Send the classification information to the RANK == 0
void DNAStrandHelper::MPI_sendClusters(){

	MPI_Send(this->partitionedCluster, this->partitionedDataSize, MPI_INT, 0, 0,MPI_COMM_WORLD);
}


// Get the classification of the data points from everyone
void DNAStrandHelper::MPI_recieveClusters(int *cluster){
	int *partitioned = this->partitionCluster();

	int clusterIndex = 0;
	for(int i =0; i < this->myWorld; i++){
		if(i == 0){

			for(int j = 0; j < this->partitionedDataSize; j++){
				cluster[clusterIndex] = this->partitionedCluster[j];
				clusterIndex++;
			}
		} else {

			int *tempCluster = new int[partitioned[i]];
			MPI_Recv(tempCluster, partitioned[i], MPI_INT, i, 0, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
			for(int j = 0; j < partitioned[i]; j++){
				cluster[clusterIndex] = tempCluster[j];
				clusterIndex++;
			}
			delete[] tempCluster;
		}


	}
}

// The slaves will recieve the partitioned data
void DNAStrandHelper::MPI_recievePartitionData(){

	MPI_Recv(&this->partitionedDataSize, 1, MPI_INT, 0, 0, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
	this->partitionedData = new char*[this->partitionedDataSize];

	for(int j = 0; j < this->partitionedDataSize; j++){
		this->partitionedData[j] = new char[this->length];
		MPI_Recv(this->partitionedData[j], this->length, MPI_CHAR, 0, 0, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
	}


}


// We are partitioning the total data in the clusters
// So that we could use this information later on when we are sending the data
int* DNAStrandHelper::partitionCluster(){
	int *partition= new int[this->myWorld];
	int overall = 0;
	for(int i = 0; i < this->myWorld; i++){
		if(i != (this->myWorld-1)){
			partition[i] = this->dataSize/this->myWorld;
			overall = overall + this->dataSize/this->myWorld;
		} else {
			partition[i] = this->dataSize - overall;
		}

	}
	return partition;
}

// Return its rank.
int DNAStrandHelper::MPI_getRank(){
	return this->myRank;
}
