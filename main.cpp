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

//int binarySearch(int arr[], int l, int r, int x) {
//    if (r >= l) {
//        int mid = (l + r - 1) / 2;
//
//        if (arr[mid] == x){
//            return mid;
//        }
//
//        if (arr[mid] < x) {
//            return binarySearch(arr, mid + 1, r, x);
//        }else {
//            return binarySearch(arr, l, mid - 1, x);
//        }
//    }
//
//
//    return -1;
//}



bool compare(int i, int j) {return (i<j);}

int** generateData(int numArray, int maxArrayLen, int maxValue, int* eachLen){
    // generate sorted data
    
    srand((unsigned) time(NULL));
    
    int** data = 0;
    data = new int* [numArray];
    
    for (int i = 0; i < numArray; i++){
        int arrayLen = rand() % (maxArrayLen) + 1;
        eachLen[i] = arrayLen;
        
        data[i] = new int[maxArrayLen];
        for (int j = 0; j < arrayLen; j++){
            data[i][j] = rand() % (maxValue) + 1;
        }
            
        if (arrayLen > 1) std::sort(data[i], data[i] + arrayLen, compare);
        
        for (int j = arrayLen; j < maxArrayLen; j++){
            data[i][j] = -1;
        }
    }
    
//    for (int i = 0; i < numArray; i++){
//        printf("%d ", eachLen[i]);
//    }
    
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
        // printf("max%d len%d each%d each+%d \n", maxValue, eachLen1[i], eachLen[i], eachLen1[i+1]);
    }

    newLen = eachLen1;
    
//    for (int i = 0; i < numArray; i++){
//        printf("newLen%d ", newLen[i]);
//    }
//    printf("\n");
    
    return maxValue;
}


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
//        printf("i:%d\n", i);
//        for (int j = 0; j < maxLen; j++){
//            printf("%d[%d, %d]\t", items[i][j].value, items[i][j].thisPos, items[i][j].nextPos);
//        }
//        printf("\n");
    }
    
    return items;
}


void fractionCascading(int** data, int numArray, int* eachLen, int** mergedData, int** thisPosData, int** nextPosData){
    int maxLen = 0;
    for (int i = 0; i < numArray; i++){
        maxLen += eachLen[i] / 2;
    }
    
    // update the last row
    int lastLen = eachLen[numArray - 1];
    mergedData[numArray - 1] = new int[maxLen];
    thisPosData[numArray - 1] = new int[maxLen];
    nextPosData[numArray - 1] = new int[maxLen];
    
    for (int j = 0; j < maxLen; j++){
        if (j < lastLen) {
            mergedData[numArray - 1][j] = data[numArray - 1][j];
            thisPosData[numArray - 1][j] = j;
        }else{
            mergedData[numArray - 1][j] = -1;
            thisPosData[numArray - 1][j] = 0;
        }
        nextPosData[numArray - 1][j] = 0;
    }
    
    // update previous rows
    for (int i = numArray - 2; i >= 0; i++){
        mergedData[i] = new int[maxLen];
        thisPosData[i] = new int[maxLen];
        nextPosData[i] = new int[maxLen];
    }
    

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
    // printf("idx: %d, thisPos: %d, nextPos: %d\n", idx, items[0][idx].thisPos, nextIdx);
    
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
    
    int key = data[0][2] + 1;
    
    int N = 1000000000;
    int* idx;
    std::clock_t start;
    double duration;

    start = std::clock();
    int i = 0;
    for (i = 0; i < N; i++){
        // if (( std::clock() - start ) / (double) CLOCKS_PER_SEC == 1.0) {break;}
        key = rand() % 100 + 1;
        idx = findPos(items, data, numArray, newLen, eachLen, key);
    }
    duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;

    std::cout<<"printf: "<< duration * 1000 << " ms'\n";
    
    printf("\nkey is: %d\n", key);
    printf("Values:\n");
    for (int i = 0; i < numArray; i++){
        printf("%d ", idx[i]);
    }
    
    return 0;
}
