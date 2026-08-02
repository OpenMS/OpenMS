/*
 *   Copyright (C) 2015-2020 Mateusz Łącki and Michał Startek.
 *
 *   This file is part of IsoSpec.
 *
 *   IsoSpec is free software: you can redistribute it and/or modify
 *   it under the terms of the Simplified ("2-clause") BSD licence.
 *
 *   IsoSpec is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 *   You should have received a copy of the Simplified BSD Licence
 *   along with IsoSpec.  If not, see <https://opensource.org/licenses/BSD-2-Clause>.
 */

#pragma once

#include <iostream>
#include <vector>
#include <cstring>
#include "isoMath.h"

namespace IsoSpec
{

inline double combinedSum(
    const int* conf, const std::vector<double>** valuesContainer, int dimNumber
){
    double res = 0.0;
    for(int i = 0; i < dimNumber; i++)
        res += (*(valuesContainer[i]))[conf[i]];
    return res;
}

inline int* getConf(void* conf)
{
    return reinterpret_cast<int*>(
        reinterpret_cast<char*>(conf) + sizeof(double)
    );
}

inline double getLProb(void* conf)
{
    double ret = *reinterpret_cast<double*>(conf);
    return ret;
}


inline double unnormalized_logProb(const int* conf, const double* logProbs, int dim)
{
    double  res = 0.0;

    for(int i = 0; i < dim; i++)
        res += minuslogFactorial(conf[i]) + conf[i] * logProbs[i];

    return res;
}

inline double calc_mass(const int* conf, const double* masses, int dim)
{
    double res = 0.0;

    for(int i = 0; i < dim; i++)
    {
        res += conf[i] * masses[i];
    }

    return res;
}



template<typename T> void printArray(const T* array, int size, const char* prefix = "")
{
    if (strlen(prefix) > 0)
        std::cout << prefix << " ";
    for (int i = 0; i < size; i++)
        std::cout << array[i] << " ";
    std::cout << std::endl;
}

template<typename T> void printVector(const std::vector<T>& vec)
{
    printArray<T>(vec.data(), vec.size());
}

template<typename T> void printOffsets(const T** array, int size, const T* offset, const char* prefix = "")
{
    if (strlen(prefix) > 0)
        std::cout << prefix << " ";
    for (int i = 0; i < size; i++)
        std::cout << array[i] - offset << " ";
    std::cout << std::endl;
}

template<typename T> void printNestedArray(const T** array, const int* shape, int size)
{
    for (int i = 0; i < size; i++)
        printArray(array[i], shape[i]);
    std::cout << std::endl;
}

//! Quickly select the n'th positional statistic, including the weights.
void* quickselect(const void** array, int n, int start, int end);


template <typename T> inline static T* array_copy(const T* A, int size)
{
    T* ret = new T[size];
    memcpy(ret, A, size*sizeof(T));
    return ret;
}

template <typename T> static T* array_copy_nptr(const T* A, int size)
{
    if(A == nullptr)
        return nullptr;
    return array_copy(A, size);
}

template<typename T> void dealloc_table(T* tbl, int dim)
{
    for(int i = 0; i < dim; i++)
    {
        delete tbl[i];
    }
    delete[] tbl;
}

template<typename T> void realloc_append(T** array, T what, size_t old_array_size)
{
    T* newT = new T[old_array_size+1];
    memcpy(newT, *array, old_array_size*sizeof(T));
    newT[old_array_size] = what;
    delete[] *array;
    *array = newT;
}

}  // namespace IsoSpec
