//
//  main.cpp
//  FractionalCascading
//
//  Created by Vincent on 2020/3/2.
//  Copyright © 2020 Vincent. All rights reserved.
//

#include <iostream>
#include <array>
#include <fstream>
#include <vector>
#include <sstream>
#include <algorithm>
#include <iterator>


#include <cstdio>
#include <ctime>

using namespace std;


bool compare(int i, int j) {return (i<j);}

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


void mergeData(int** mergedData){
    
}

typedef struct {
    int value = -1;
    int thisPos = 0;
    int nextPos = 0;
}item;

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

int binarySearchArr(int arr[], int arraySize, int key) {

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

    }

    newLen = eachLen1;
    
    return maxValue;
}


int findIdxInArray(int key, int size, int* array){
    
    int pivot, left = 0, right = size - 1;
    while (left <= right) {
      pivot = left + (right - left) / 2;
      if (array[pivot] == key) return pivot;
      if (key < array[pivot]) 
      right = pivot - 1;
      else left = pivot + 1;
    }
    return left;
}

int findIdxInItem(int key, int size, item* items){
    int pivot, left = 0, right = size - 1;
    while (left <= right) {
      pivot = left + (right - left) / 2;
      if (items[pivot].value == key) 
      return pivot;
      if (key < items[pivot].value) 
      right = pivot - 1;
      else left = pivot + 1;
    }
    return left;
}


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

int* findPos(item** items, int** data, int numArray, int* newLen, int* oriLen, int key){

    int* valueArray = new int[numArray];
    
    // find idx of first row
    int idx = binarySearchItem(items[0], newLen[0], key);
    if (idx > newLen[0] - 1){
        idx = newLen[0] - 1;
    }
    
    int thisIdx = items[0][idx].thisPos;
    if (thisIdx > oriLen[0] - 1){
        thisIdx = oriLen[0] - 1;
    }
    
    valueArray[0] = data[0][thisIdx];

    int nextIdx = items[0][idx].nextPos;
    if (nextIdx > newLen[1] - 1) {
        nextIdx = newLen[1] - 1;
    }
    
    // other lines
    for (int i = 1; i < numArray; i++){
        item selectedItem = items[i][nextIdx];
        
        if ((nextIdx > 0) && items[i][nextIdx - 1].value > key){
            selectedItem = items[i][nextIdx - 1];
        }
        
        if ((nextIdx < newLen[i] - 1) && items[i][nextIdx + 1].value < key){
            selectedItem = items[i][nextIdx + 1];
        }
        
        thisIdx = selectedItem.thisPos;
        if (thisIdx > oriLen[i] - 1){
            thisIdx = oriLen[i] - 1;
        }
        
        valueArray[i] = data[i][thisIdx];
        
        nextIdx = selectedItem.nextPos;
        if (i < numArray - 1 && nextIdx > newLen[i + 1] - 1){
            nextIdx = newLen[i + 1] - 1;
        }
        
    }
    
    return valueArray;
}

int main(int argc, const char * argv[]) {


    
    int numArray = 10;
    int maxArrayLen = 5;
    int maxValue = 100;
    
    int eachLen[numArray];
    
    int** data = generateData(numArray, maxArrayLen, maxValue, eachLen);
    show2dData(data, numArray, maxArrayLen);
    
    item** items = 0;
    int* maxLength = new int[1];
    int* newLen = new int[numArray];
    items = fractionalCascading(data, numArray, eachLen, maxLength, newLen);
    
    // for loop in fc function not get in
    printf("Data after fractional cascading:\n");
    for (int i = 0; i < numArray; i++){
        for (int j = 0; j < *maxLength; j++){
            printf("%3d[%d, %d]\t\t", items[i][j].value, items[i][j].thisPos, items[i][j].nextPos);
        }
        printf("\n");
    }
    
    int* idx;
    int key;
    
    std::clock_t start;
    double duration;

    start = std::clock();
	
    // change number of iterations here
	int N = 100000;
    for(int i = 0; i< N; i++){
	  key = rand() % 100 + 1;
      idx = findPos(items, data, numArray, newLen, eachLen, key);
    }
      
	duration = ( (std::clock() - start) * 1.0 ) / (double) CLOCKS_PER_SEC * 1000;
	printf("Time %f\n", duration);
       

    
    printf("\nkey is: %d\n", key);
    printf("Values:\n");
    for (int i = 0; i < numArray; i++){
        printf("%d ", idx[i]);
    }
    
    return 0;
}


