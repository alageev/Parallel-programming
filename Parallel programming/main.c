//
//  main.cpp
//  Parallel programming
//
//  Created by Алексей Агеев on 12.02.2021.
//

#include <math.h>
#include <stdio.h>
#include <sys/time.h>
#include <stdlib.h>


//MARK:- Task
// A = сount(Агеев) * сount(Алексей) * сount(Дмитриевич) = 5 * 7 * 10 = 350
int const aTask = 350;

// X_1 = 1 + ((350 mod 47) mod 7) = 1 + (21 mod 7) = 1
// Гиперболический синус с последующим возведением в квадрат

// X_2 = 1 + ((350 mod 47) mod 8) = 1 + (21 mod 8) = 6
// Десятичный логарифм, возведенный в степень e

// X_3 = 1 + ((350 mod 47) mod 6) = 1 + (21 mod 6) = 4
// Выбор большего (т.е. M2[i] = max(M1[i],M2[i])))

// X_4 = 1 + ((350 mod 47) mod 7) = 1 + (21 mod 7) = 1
//Сортировка выбором (Selection sort)

void generateTwoArrays(double*, double*, int, unsigned int*);
void map(double*, double*, int);
void merge(double*, double*, int);
void selectionSort(double*, const int);
void reduce(double*, double*, const int);


//MARK:- main()
int main(int argc, const char * argv[]) {
    
    int const length = atoi(argv[1]);
    unsigned seed;
    double results[50] = { 0 };

    struct timeval startTime;
    struct timeval endTime;
    gettimeofday(&startTime, NULL);

    for (unsigned i = 0; i < 50; i++) {
        srand(i);
        double firstArray[length];
        double secondArray[length / 2];

        //MARK:- Generate
        seed = i;
        generateTwoArrays(firstArray, secondArray, length, &seed);


        //MARK:- Map
        map(firstArray, secondArray, length);


        //MARK:- Merge
        merge(firstArray, secondArray, length);


        //MARK:- Sort
        selectionSort(secondArray, length / 2);


        //MARK:- Reduce
        reduce(secondArray, results + i, length / 2);
    }

    gettimeofday(&endTime, NULL);
    long const workTime = 1000 * (endTime.tv_sec - startTime.tv_sec) + (endTime.tv_usec - startTime.tv_usec) / 1000;

    printf("\nN=%d. Milliseconds passed: %ld\n", length, workTime);
}

/// This function fills two arrays using generation task
/// @param first First array
/// @param second Second Array
/// @param length Length of the first array
/// @param seed the seed for rand_r() function
void generateTwoArrays(double* first, double* second, int length, unsigned* seed) {
    #pragma omp parallel
    {
        #pragma omp for
        for (int j = 0; j < length; j++) {
            first[j] = rand_r(seed) % aTask + 1;
        }
        #pragma omp for
        for (int j = 0; j < length / 2; j++) {
            second[j] = aTask + rand_r(seed) % (9 * aTask + 1);
        }
    }
}

/// This function maps some functions over the array
/// @param first First array
/// @param second Second Array
/// @param length Length of the first array
void map(double* first, double* second, int length) {
    
    #pragma omp parallel
    {
        #pragma omp for
        for (int j = 0; j < length; j++) {
            first[j] = pow(sinh(first[j]), 2);
        }
        
        #pragma omp for//??????????????
        for (int j = length / 2 - 1; j > 0; j--) {
            second[j] = pow(log10(second[j] + second[j - 1]), M_E);
        }
    }
}

/// This function merges two arrays
/// @param first First array
/// @param second Second Array
/// @param length Length of the first array
void merge(double* first, double* second, int length) {
    
    #pragma omp parallel for shared(first, second)
    for (int j = 0; j < length / 2; j++) {
        if (first[j] > second[j]) {
            second[j] = first[j];
        }
    }
}

/// This function sorts given array using selection sorting algorithm
/// @param array array that should be sorted
/// @param length length of the array
void selectionSort(double* array, const int length) {
    for (int i = 0; i < length - 1; i++) {
        int minIndex = i;
        
        #pragma omp parallel for shared(minIndex)
        for (int j = i + 1; j < length; j++) {
            if (array[j] < array[minIndex]) {
                #pragma omp atomic
                minIndex = j;
            }
        }

        if (minIndex != i) {
            double const temp = array[i];
            array[i] = array[minIndex];
            array[minIndex] = temp;
        }
    }
}

/// This function
/// @param array array that shoud be reduced
/// @param result a value that should be returned as result of this function
/// @param length length of the array
void reduce(double* array, double* result, const int length) {

    double minNonZero = __DBL_MAX__;
    double res = 0;
    
    #pragma omp parallel
    {
        #pragma omp for reduction(min:minNonZero)
        for (int j = 0; j < length / 2; j++) {
            double value = array[j];
            if (minNonZero > value && value != 0) {
                minNonZero = value;
            }
        }
        
        

        #pragma omp for reduction(+:res)
        for (int j = 0; j < length / 2; j++) {
            if ((int)floor(array[j] / minNonZero) % 2) {
                res += sin(array[j]);
            }
        }
    }
    *result = res;
}
