//
//  main.cpp
//  FractionalCascading
//
//  Created by Vincent on 2020/3/2.
//  Copyright © 2020 Vincent. All rights reserved.
//

#include <array>
#include <fstream>
#include <vector>
#include <sstream>
#include <algorithm>
#include <iterator>
#include <curand.h>
#include "curand_kernel.h"
#include <time.h>

using namespace std;

bool compare(int i, int j) {return (i<j);}

__host__
vector<int> GenerateDiffNumber(int min,int max,int num)
{
    int rnd;
    vector<int> diff;
    vector<int> tmp;

    for(int i = min;i < max+1 ; i++ )
    {
        tmp.push_back(i);
    }
    srand((unsigned)time(0));
    for(int i = 0 ; i < num ; i++)
    {
        do{
            rnd = min+rand()%(max-min+1);
     
        }while(tmp.at(rnd-min)==-1);
        diff.push_back(rnd);
        tmp.at(rnd-min) = -1;
    }
    return diff;
}

__host__
int** generateData(int row, int col, int maxValue, int* eachLen){
    // generate sorted data
    
    srand((unsigned) time(NULL));
    
    int** data = 0;
    data = new int* [row];
    
    vector<int> totalRandons = GenerateDiffNumber(1, 500, 200);
    int k = 0;

    for (int i = 0; i < row; i++){
        int actualLen = rand() % (col) + 1;//length of nonzero elements in each row
        eachLen[i] = actualLen;
        
        data[i] = new int[col];

        
         for (int j = 0; j < actualLen; j++){
           // data[i][j] = rand() % (maxValue) + 1;
           data[i][j] = totalRandons[k++];
        }
            
        if (actualLen > 1) std::sort(data[i], data[i] + actualLen, compare);
        
        for (int j = actualLen; j < col; j++){
            data[i][j] = -1;
        }
    }
    
    return data;
}

typedef struct {
    int value = -1;
    int thisPos = 0;
    int nextPos = 0;
}item;

__device__
void binarySearchGPU(int* arr, int arraySize, int key, int* &pos) {

	printf("here bsn\n");

	int l = 0;
    int r = arraySize - 1;
    int mid = (l+r)/2;
    
    while(l <= r){
        mid = (l + r)/2;
        if(arr[mid] == key){
            *pos = mid;
			break;
        }
        
        if(arr[mid] > key){
            r = mid - 1;
        }else {
            l = mid + 1;
        }
    }
    
    *pos = l;
}

__host__ __device__
int binarySearchItem(item* items, int arraySize, int key){
    
    int l = 0;
    int r = arraySize - 1;
    int mid = (l+r)/2;
    
    while(l <= r){
        mid = (l + r)/2;
        if(items[mid].value == key){
            return mid;
        }
        
        if(items[mid].value > key){
            r = mid - 1;
        }else {
            l = mid + 1;
        }
    }
    
    return l;
}


__host__ __device__
int binarySearchArr(int* arr, int arraySize, int key) {

    int l = 0;
    int r = arraySize - 1;

    while(l <= r) {
        int mid = (l+r)/2;
        if (arr[mid] == key){
            return mid;
        }

        if (arr[mid] > key) {
            r = mid - 1;
        }else {
            l = mid + 1;
        }
    }

    return l;
}

__host__
int findMaxLen(int numArray, int* eachLen, int* &newLen){

    if (numArray == 0) return 0;
    
    int maxValue = 0;

    int* eachLen1 = new int[numArray];
    for (int i = 0; i < numArray; i++){
        eachLen1[i] = eachLen[i];
    }

    // calculate the maxLen: iterate from bottom back to top
    // this len = thisLen + preLen/2
    for (int i = numArray - 2; i >= 0; i--){
        eachLen1[i] = eachLen[i] + eachLen1[i+1] / 2;
        if (eachLen1[i] > maxValue) maxValue = eachLen1[i];
        // printf("max%d len%d each%d each+%d \n", maxValue, eachLen1[i], eachLen[i], eachLen1[i+1]);
    }

    newLen = eachLen1;
    
//    for (int i = 0; i < numArray; i++){
//        printf("newLen%d ", newLen[i]);
//    }
//    printf("\n");
    
    return maxValue;
}

__host__
int findIdxInArray(int key, int size, int* array){
    
    int pivot, left = 0, right = size - 1;
    while (left <= right) {
      pivot = left + (right - left) / 2;
      if (array[pivot] == key) return pivot;
      if (key < array[pivot]) right = pivot - 1;
      else left = pivot + 1;
    }
    return left;
}

__host__
int findIdxInItem(int key, int size, item* items){
    int pivot, left = 0, right = size - 1;
    while (left <= right) {
      pivot = left + (right - left) / 2;
      if (items[pivot].value == key) return pivot;
      if (key < items[pivot].value) right = pivot - 1;
      else left = pivot + 1;
    }
    return left;
}

__host__
item** fractionalCascading(int** data, int numArray, int* eachLen, int* maxLength, int* &newLen){
    
    int maxLen = findMaxLen(numArray, eachLen, newLen);
    maxLength[0] = maxLen;
    item** items = 0;
    items = new item* [numArray];
    
    // initialize the last row
    int lastLen = eachLen[numArray - 1];
    items[numArray - 1] = new item[maxLen];
    for (int i = 0; i < maxLen; i++){
        if (i < lastLen){
            items[numArray - 1][i].value = data[numArray - 1][i];
            items[numArray - 1][i].thisPos = i;
            items[numArray - 1][i].nextPos = 0;
        }
    }
    

    // other rows, merge the next with this one
    for (int i = (numArray - 2); i >= 0; i--){
        
        // this line i, next line i + 1
        
        items[i] = new item[maxLen];
        int thisLen = eachLen[i];
        int nextLen = newLen[i + 1];
        int thisNewLen = thisLen + nextLen / 2;
        
        // x for this line, y for next line, z for new item array
        // place element first, then find positions
        int x = 0, y = 1, z = 0;
        while(x < thisLen && y < nextLen){
            if (z < thisNewLen) {
                if (data[i][x] <= items[i + 1][y].value){
                    items[i][z].value = data[i][x];
                    items[i][z].thisPos = x;
                    items[i][z].nextPos = binarySearchItem(items[i + 1], newLen[i + 1], data[i][x]);
                    x++;
                    z++;
                }else {
                    items[i][z].value = items[i + 1][y].value;
                    items[i][z].thisPos = binarySearchArr(data[i], thisLen, items[i + 1][y].value);
                    items[i][z].nextPos = y;
                    y += 2;
                    z++;
                }
            }
        }
        
        while (x < thisLen) {
            items[i][z].value = data[i][x];
            items[i][z].thisPos = x;
            items[i][z].nextPos = findIdxInItem(data[i][x], nextLen, items[i+1]);
            x++;
            z++;
        }
        
        while (y < nextLen) {
            items[i][z].value = items[i + 1][y].value;
            items[i][z].thisPos = findIdxInArray(items[i + 1][y].value, thisLen, data[i]);
            items[i][z].nextPos = y;
            y += 2;
            z++;
        }
    }
    
    return items;
}

__host__
void show2dData(int** data, int numArray, int maxArrayLen){
    printf("Sorted Raw Data:\n");
    for (int i = 0; i < numArray; i++){
        for (int j = 0; j < maxArrayLen; j++){
            printf("%d\t", data[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}


const int numArray = 10;
const int maxArrayLen = 5;
const int maxValue = 100;
const int length = 15;


__global__
void kernelSetRandom(curandState *curandStates, int N, long clock_for_rand){
	int idx = threadIdx.x + blockIdx.x * blockDim.x;
	if (idx > 1000) return;
	// int stride = blockDim.x * gridDim.x;
	
	// if (idx > num)	return;
	for (int a = 0; a < 1000; a ++){
		curand_init(clock_for_rand, idx, 0, &curandStates[a]);
	}
}

__global__
void findPos(int** valueArr, int** thisPosArr, int** nextPosArr, int** data, int numArray, int* newLen, int* oriLen, int* valueArray, int N, int* keyOut, curandState* curand_states){
    
    // find idx of first row
	// int idx = binarySearchArr(valueArr[0], newLen[0], key);
	/*
	if (threadIdx.x==0){
		for(int i = 0; i < numArray; i++){
			for(int j = 0; j < length; j++){
				printf("%d ", valueArr[i][j]);
			}
		}
	}
	*/

	int threadIndex = threadIdx.x + blockIdx.x * blockDim.x;
	int stride = blockDim.x * gridDim.x;

	for (int a = threadIndex; a < N; a += stride){
		int b = a % 1000;
	// if (threadIndex > N)	return;
	curandState localState = curand_states[b];
	int key = (int)abs(curand(&localState)) % maxValue + 1;
	// printf("inside key: %d\n", key);
	*keyOut = key;
	
	int idx;
    int l = 0;
    int r = newLen[0] - 1;

    while(l <= r) {
        int mid = (l+r)/2;
        if (valueArr[0][mid] == key){
            idx = mid;
        }

        if (valueArr[0][mid] > key) {
            r = mid - 1;
        }else {
            l = mid + 1; 
        }
    }

    idx = l;
	// end of binary search	
	
    if (idx > newLen[0] - 1){
        idx = newLen[0] - 1;
    }
    
    int thisIdx = thisPosArr[0][idx];
    if (thisIdx > oriLen[0] - 1){
        thisIdx = oriLen[0] - 1;
    }
    
	valueArray[0] = data[0][thisIdx];
	//printf("thisIdx %d ", thisIdx); 

    int nextIdx = nextPosArr[0][idx];
    if (nextIdx > newLen[1] - 1) {
        nextIdx = newLen[1] - 1;
    }
	 

    // other lines
    for (int i = 1; i < numArray; i++){

		int selectedItemThisPos = thisPosArr[i][nextIdx];
		int selectedItemNextPos = nextPosArr[i][nextIdx];
		// printf("i %d, nextIdx %d, value %d\n", i, nextIdx, selectedItemValue);
        
        if ((nextIdx > 0) && valueArr[i][nextIdx - 1] > key){
			selectedItemThisPos = thisPosArr[i][nextIdx - 1];
			selectedItemNextPos = nextPosArr[i][nextIdx - 1];
        }
        
        if ((nextIdx < newLen[i] - 1) && valueArr[i][nextIdx + 1] < key){
			selectedItemThisPos = thisPosArr[i][nextIdx + 1];
			selectedItemNextPos = nextPosArr[i][nextIdx + 1];
        }
        
        thisIdx = selectedItemThisPos;
        if (thisIdx > oriLen[i] - 1){
            thisIdx = oriLen[i] - 1;
        }
        
        valueArray[i] = data[i][thisIdx];
        
        nextIdx = selectedItemNextPos;
        if ((i < numArray - 1) && (nextIdx > newLen[i + 1] - 1)){
            nextIdx = newLen[i + 1] - 1;
        }   
    }
}


/*
	for (int i = 0; i < numArray; i++){
		printf("%d ", valueArray[i]);
	}
*/
}

int main(int argc, const char * argv[]) {
    // eachLen (oriLen)    
    int* eachLen = new int[numArray];
	int* eachLenK;
	cudaMallocManaged(&eachLenK, numArray * sizeof(int));
    
	// data
    int** data = generateData(numArray, maxArrayLen, maxValue, eachLen);
	cudaMemcpy(eachLenK, eachLen, numArray * sizeof(int), cudaMemcpyHostToDevice);
    show2dData(data, numArray, maxArrayLen);

	int dataSize = numArray*maxArrayLen;
	int* dataK; cudaMallocManaged(&dataK, dataSize * sizeof(int));
	std::copy(&data[0][0], &data[0][0] + dataSize, dataK);

	// new length
    int* newLen = new int[numArray];
	int* newLenK;
	cudaMallocManaged(&newLenK, numArray * sizeof(int));

	// items
	item** items = 0;
    int* maxLength = new int[1];
    items = fractionalCascading(data, numArray, eachLen, maxLength, newLen);
	cudaMemcpy(newLenK, newLen, numArray * sizeof(int), cudaMemcpyHostToDevice);
    
    // print items on screen
    printf("Data after fractional cascading:\n");
    for (int i = 0; i < numArray; i++){
        for (int j = 0; j < *maxLength; j++){
            printf("%3d[%d, %d] ", items[i][j].value, items[i][j].thisPos, items[i][j].nextPos);
        }
        printf("\n");
    }
    
	// get key
    // int key = data[0][2] + 1;
	
	// convert items struct to arrays
	// const int length = *maxLength;

	int** valueArr = new int* [numArray];
	int** thisPosArr = new int* [numArray];
	int** nextPosArr = new int* [numArray];
	for (int i = 0; i < numArray; i++) {
		valueArr[i] = new int[length];
		thisPosArr[i] = new int[length];
		nextPosArr[i] = new int[length];
		for (int j = 0; j < length; j++) {
			if (j < *maxLength) {
				valueArr[i][j] = items[i][j].value;
				thisPosArr[i][j] = items[i][j].thisPos;
				nextPosArr[i][j] = items[i][j].nextPos;
			}else{
				valueArr[i][j] = -1;
				thisPosArr[i][j] = 0;
				nextPosArr[i][j] = 0;
			}
		}	
	}



	// convert items arrays into cuda format
	// testing for new mem allocation
	int* data_valueArr = (int *)malloc(sizeof(int) * length * numArray);
	int* data_thisPosArr = (int *)malloc(sizeof(int) * length * numArray);
	int* data_nextPosArr = (int *)malloc(sizeof(int) * length * numArray);

	for (int i = 0; i < numArray; i++){
		for(int j = 0; j < length; j++){
		data_valueArr[i * length + j] = valueArr[i][j];
		data_thisPosArr[i * length + j] = thisPosArr[i][j];
		data_nextPosArr[i * length + j] = nextPosArr[i][j]; 
		}
	}
	
	int* data_data = (int *)malloc(sizeof(int) * length * numArray);
	for (int i = 0; i < numArray; i++) {
		for (int j = 0; j < maxArrayLen; j++){
			data_data[i * maxArrayLen + j] = data[i][j];
		}
	}


	int** d_valueArr;
	int* d_data_valueArr;
	int** d_thisPosArr;
	int* d_data_thisPosArr;
	int** d_nextPosArr;
	int* d_data_nextPosArr;
	int** d_data;
	int* d_data_data;
	cudaMalloc((void**)&d_valueArr, sizeof(int *) * numArray);
	cudaMalloc((void**)&d_data_valueArr, sizeof(int) * numArray * length);

	cudaMalloc((void**)&d_thisPosArr, sizeof(int *) * numArray);
	cudaMalloc((void**)&d_data_thisPosArr, sizeof(int) * numArray * length);

	cudaMalloc((void**)&d_nextPosArr, sizeof(int *) * numArray);
	cudaMalloc((void**)&d_data_nextPosArr, sizeof(int) * numArray * length);

	cudaMalloc((void**)&d_data, sizeof(int *) * numArray);
	cudaMalloc((void**)&d_data_data, sizeof(int) * numArray * maxArrayLen);

	for (int i = 0; i < numArray; i++){
		valueArr[i] = d_data_valueArr + length * i;
		thisPosArr[i] = d_data_thisPosArr + length * i;
		nextPosArr[i] = d_data_nextPosArr + length * i;
		data[i] = d_data_data + maxArrayLen * i;
	}

	cudaMemcpy(d_valueArr, valueArr, sizeof(int*) * numArray, cudaMemcpyHostToDevice);
	cudaMemcpy(d_data_valueArr, data_valueArr, sizeof(int) * numArray * length, cudaMemcpyHostToDevice);
	cudaMemcpy(d_thisPosArr, thisPosArr, sizeof(int*) * numArray, cudaMemcpyHostToDevice);
	cudaMemcpy(d_data_thisPosArr, data_thisPosArr, sizeof(int) * numArray * length, cudaMemcpyHostToDevice);
	cudaMemcpy(d_nextPosArr, nextPosArr, sizeof(int*) * numArray, cudaMemcpyHostToDevice);
	cudaMemcpy(d_data_nextPosArr, data_nextPosArr, sizeof(int) * numArray * length, cudaMemcpyHostToDevice);
	cudaMemcpy(d_data, data, sizeof(int*) * numArray, cudaMemcpyHostToDevice);
	cudaMemcpy(d_data_data, data_data, sizeof(int) * numArray * maxArrayLen, cudaMemcpyHostToDevice);

	int N = 1000000;
	int blockSize = 1024;
    int numBlock = (N + blockSize - 1) / blockSize;

	int* resultK;
	cudaMallocManaged(&resultK, numArray * sizeof(int));
	
	// prepare for kernel random number generation (key for fc)
	curandState* dev_states;
	cudaMalloc((void**) &dev_states, sizeof(curandState) * 1000);
	long clock_for_rand = clock();
	kernelSetRandom<<<numBlock, blockSize>>>(dev_states, N, clock_for_rand);

	// kernel FC function
	int* d_key;
	cudaMallocManaged(&d_key, sizeof(int));

	cudaEvent_t start, stop;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);
	float elapsedTime;
	cudaEventRecord(start, 0);
	
    findPos<<<numBlock, blockSize>>>(d_valueArr, 
									 d_thisPosArr, 
								     d_nextPosArr, 
									 d_data, 
									 numArray, newLenK, eachLenK, resultK, N, d_key, dev_states);
  	

	cudaEventRecord(stop, 0);
	cudaEventSynchronize(stop);
	cudaEventElapsedTime(&elapsedTime, start, stop);
	printf("Time: %f ms\n", elapsedTime);
	cudaEventDestroy(start);
	cudaEventDestroy(stop);
	
	cudaDeviceSynchronize();

    
	int* h_key;
	h_key = (int *)malloc(sizeof(int));
	cudaMemcpy(h_key, d_key, sizeof(int), cudaMemcpyDeviceToHost);
    printf("\nkey: %d\n", h_key[0]);
	
	int* hostResult = new int[numArray];
	hostResult = (int *)malloc(numArray * sizeof(int));
	cudaMemcpy(hostResult, resultK, numArray * sizeof(int), cudaMemcpyDeviceToHost);
    printf("Values: ");
    for (int i = 0; i < numArray; i++){
        printf("%d ", hostResult[i]);
    }

	printf("\n"); 

	cudaFree(eachLen);
	cudaFree(newLenK);
	cudaFree(d_valueArr);
	cudaFree(d_data_valueArr);
	cudaFree(d_thisPosArr);
	cudaFree(d_data_thisPosArr);
	cudaFree(d_nextPosArr);
	cudaFree(d_data_nextPosArr);
	cudaFree(d_data);
	cudaFree(d_data_data);

    
    return 0;
}
