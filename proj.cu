#include <iostream>
#include <math.h>
#include <string>


__host__
/*
int ** genData(int numArray, int maxArrayLen, int maxValue) {
	
 	int data[numArray][maxArrayLen];
	
	srand((unsigned) time (NULL));
	
	for (int i = 0; i < numArray; i++) {
		int thisArr[maxArrayLen] = [0];
		for (int j = 0; j < maxArrayLen; j++) {
			
		}	
	}
}
*/

__host__
int** readData(std::string filePath){
	std::ifstream myfile(filePath);

	if (!myfile) {
		printf("file not found");
		return NULL;
	}

	return NULL;
	


}


int main() {
	int numArray = 10;
	int maxArrayLen = 20;
	int maxValue = 500;

	int** data = 


	return 0;
}
