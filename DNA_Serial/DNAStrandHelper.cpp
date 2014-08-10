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
#include <fstream>
#include <string>
using namespace std;

// Constructor
DNAStrandHelper::DNAStrandHelper(int length, int dataSize, int clusters) {
	// TODO Auto-generated constructor stub

	srand(time(0));

	this->length = length;
	this->dataSize = dataSize;
	this->clusters = clusters;

}

// Destructor
DNAStrandHelper::~DNAStrandHelper() {
}

// Function for generating a DNA Strand
char* DNAStrandHelper::getRamdomDNAStrand(){

	char* strand = new char[this->length];
	for(int i = 0 ; i < length; i++){
		strand[i] = this->getRandomType(); //strandExample[rand()%50];
	}

	return strand;
}

// Function for selecting one sequence
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

// Helper function to print out the Strand
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
			cout<< value; break;
		}
	}

	cout<<endl;
}

// Function comparing the two DNA Strand and return the distance
int DNAStrandHelper::compareTwo(char* strandA, char* strandB){
	int counter = 0;

	for(int i = 0; i < this->length; i++){
		if(strandA[i] != strandB[i]){
			counter++;
		}
	}

	return counter;
}

// Function that will pick random centroid from the Data
char** DNAStrandHelper::pickCentriods(char** data){

	cout<<"Picking centroids"<<endl;
	char** centroids = new char*[this->clusters]();

	for(int i = 0; i < this->clusters; i++){
		centroids[i] = new char[this->length]();

		int index = rand()%this->dataSize;
		for(int j = 0; j < this->length; j++){
			centroids[i][j] = data[index][j];
		}
	}
	return centroids;
}

// Function to actually do the clustering
void DNAStrandHelper::findCluster(char **data, char **centriods, int *cluster, double threshold, int maxIteration){
	cout<<"Starting to find Cluster"<<endl;
	int iteration = 0;

	double clusterUpdated = this->dataSize;
	double comparedValue = this->dataSize * threshold;
	while(   clusterUpdated > comparedValue    && maxIteration > iteration ){
		clusterUpdated = 0;
		iteration++;


		for(int i = 0; i < this->dataSize; i++){
				// this will give a validated cluster for the data point
				int validatedCluster = this->validateStrandForCluster(data[i], centriods);
				if( cluster[i] != validatedCluster){
					cluster[i] = validatedCluster;
					clusterUpdated++;
				}

			}

		// Doing the validation of the new centroids
		this->validateCentriod(data,centriods,cluster);

		cout <<"::::::: Current Changing Count = " << clusterUpdated << " Threshold = " <<  comparedValue << endl;
	}

}

// Function that will use the centroids and classify what cluster it belongs to
int DNAStrandHelper::validateStrandForCluster(char *strand, char **centriods){
	int index = 0;
	int minDistance = this->length * 2;
	for(int i = 0; i < this->clusters; i++){
		// getting the distance of the two
		int distance = this->compareTwo(centriods[i], strand);
		if( distance < minDistance){
			index = i;
			minDistance = distance;
		}

	}
	return index;
}

// This function will recalculate the centroid
void DNAStrandHelper::validateCentriod(char **data, char **centriods, int *cluster){

	// Setup a temp object that will contain the information
	double ***middleStrand = new double**[this->clusters];
	for(int mid = 0; mid < this->clusters; mid++){
		middleStrand[mid] = new double*[this->length];
		for(int len = 0; len<this->length; len++){
			for(int c = 0; c < 4; c++){
				middleStrand[mid][len] = new double[4];
				middleStrand[mid][len][c] = 0;
			}

		}
	}

	// Keep up with the cluster count
	int *clusterCounter = new int[this->clusters];
	for(int index = 0; index < this->clusters; index++){
		clusterCounter[index] = 0;
	}

	for(int i = 0; i < this->clusters; i++){
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
	for(int i = 0; i < this->clusters; i++){
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


	// Cleaning up
	for(int mid = 0; mid < this->clusters; mid++){
		for(int len = 0; len<this->length; len++){
				delete[] middleStrand[mid][len];
			}
			delete[] middleStrand[mid];
	}


	delete[] middleStrand;
	delete[] clusterCounter;
}

// Printing all data into a file
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
