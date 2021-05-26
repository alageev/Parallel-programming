//
//  main.c
//  Parallel programming
//
//  Created by Алексей Агеев on 12.02.2021.
//

#include <math.h>
#include <stdio.h>
#include <sys/time.h>
#include <stdlib.h>
#include <unistd.h>

//#include "OpenMP sources/omp.h"

#ifdef _OPENMP

#include <omp.h>

#else

double omp_get_wtime() {
    return 0.0;
}

int omp_get_num_procs() {
    return 1;
}


#endif// ifdef _OPENMP



//MARK:- Task
// A = сount(Агеев) * сount(Алексей) * сount(Дмитриевич) = 5 * 7 * 10 = 350
int const aTask = 350;
#define NUMBER_OF_ITERATIONS 50

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
void selectionSortOfPart(double*, const int, const int);
//void mergeParts(double*, const int, const int);
void mergeParts(double*, const int, const int, const int);
void timer(double*);

//MARK:- main()
int main(int argc, const char * argv[]) {
    double completion = 0;

    #pragma omp parallel sections shared(completion)
    {
        
        #pragma omp section
        {
            timer(&completion);
        }

        #pragma omp section
        {
            int const length = atoi(argv[1]);
            unsigned seed;
            double results[NUMBER_OF_ITERATIONS] = { 0 };

            double const startTime = omp_get_wtime();

            for (unsigned i = 0; i < NUMBER_OF_ITERATIONS; i++) {
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

                completion += 100.0 / NUMBER_OF_ITERATIONS;
            }

            printf("\nN=%d. Milliseconds passed: %f\n", length, 1000 * (omp_get_wtime() - startTime));
        }
    }
}

/// This function fills two arrays using generation task
/// @param first First array
/// @param second Second Array
/// @param length Length of the first array
/// @param seed the seed for rand_r() function
void generateTwoArrays(double* first, double* second, int length, unsigned* seed) {
    //force single-threaded for correct result
    for (int j = 0; j < length; j++) {
        first[j] = rand_r(seed) % aTask + 1;
    }

    for (int j = 0; j < length / 2; j++) {
        second[j] = aTask + rand_r(seed) % (9 * aTask + 1);
    }
}

/// This function maps some functions over the array
/// @param first First array
/// @param second Second Array
/// @param length Length of the first array
void map(double* first, double* second, int length) {

    #pragma omp parallel for
    for (int j = 0; j < length; j++) {
        first[j] = pow(sinh(first[j]), 2);
    }

    for (int j = length / 2 - 1; j > 0; j--) {
        second[j] = pow(log10(second[j] + second[j - 1]), M_E);
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
    int const numberOfSortingThreads = omp_get_num_procs() - 1;
    int const partSize = length / numberOfSortingThreads;
    int const lastPartSize = length - (numberOfSortingThreads) * partSize;

    #pragma omp parallel sections
    {
        for (int i = 0; i < numberOfSortingThreads; i++) {
            #pragma omp section
            selectionSortOfPart(array, i * partSize, partSize);
        }

        #pragma omp section
        selectionSortOfPart(array, numberOfSortingThreads * partSize, lastPartSize);
    }

    mergeParts(array, partSize, length, numberOfSortingThreads + 1);
}

/// This function reduces array to result value
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


/// This function sorts given array using selection sorting algorithm
/// @param array array that should be sorted
/// @param offset offset index
/// @param length length of the array
void selectionSortOfPart(double* array, const int offset, const int length) {

    for (int i = offset; i < offset + length - 1; i++) {
        int minIndex = i;
        int localIndex = minIndex;

        #pragma omp parallel shared(minIndex) private(localIndex)
        {
            #pragma omp for
            for (int j = i + 1; j < offset + length; j++) {
                if (array[j] < array[localIndex]) {
                    localIndex = j;
                }
            }

            #pragma omp critical
            if (array[localIndex] < array[minIndex]) {
                minIndex = localIndex;
            }
        }

        if (minIndex != i) {
            double const temp = array[i];
            array[i] = array[minIndex];
            array[minIndex] = temp;
        }
    }
}


/// This function merges parts of given array in ascending order
/// @param array array whose parts should be merged
/// @param partSize size of every of 0...k-1 part
/// @param length length of the part
/// @param numberOfParts number of parts in array
void mergeParts(double* array, const int partSize, const int length, const int numberOfParts) {

    selectionSortOfPart(array, 0, length);

    // Сергей Мороз сказал, что mergesort тут использовать не надо было
    // надо было просто отсортировать массив ещё раз
    // честно плакал, когда услышал это
    // потому что дебажил этот код как никогда не дебажил
    // оставлю его тут на память

    /*
    double copy[length];
    int indices[numberOfParts];
    int indicesCeil[numberOfParts];


    #pragma omp parallel shared(copy, array)
    {
        #pragma omp for nowait
        for (int i = 0; i < length; i++) {
            copy[i] = array[i];
        }

        #pragma omp for
        for (int i = 0; i < numberOfParts; i++) {
            indices[i] = i * partSize;
        }

        #pragma omp for
        for (int i = 0; i < numberOfParts - 1; i++) {
            indicesCeil[i] = indices[i + 1];
        }

        indicesCeil[numberOfParts - 1] = length;
    }

    for (int i = 0; i < length; i++) {
        int position = 0;

        for (int j = 0; j < numberOfParts; j++) {
            if (indices[j] < indicesCeil[j]) {
                position = j;
                break;
            }
        }

        int localPosition = position;

        #pragma omp parallel shared(position) private(localPosition)
        {

            #pragma omp for
            for (int j = 0; j < numberOfParts; j++) {
                if (indices[j] < indicesCeil[j]) {
                    if (copy[indices[j]] < copy[indices[localPosition]]) {
                        localPosition = j;
                    }
                }
            }

            #pragma omp critical
            if (copy[indices[localPosition]] < copy[indices[position]]) {
                position = localPosition;
            }
        }

        array[i] = copy[indices[position]];
        indices[position]++;
    }
    */
}

void timer(double* completion) {
    int second = 0;
    while (*completion < 100) {
        sleep(1);
        second++;
        printf("%d: %f\n", second, *completion);
    }
}
