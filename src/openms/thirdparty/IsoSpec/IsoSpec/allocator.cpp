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


#include "allocator.h"

namespace IsoSpec
{

template <typename T>
Allocator<T>::Allocator(const int dim_, const int tabSize_): currentId(-1), dim(dim_), tabSize(tabSize_)
{
    currentTab = new T[dim * tabSize];
}

template <typename T>
Allocator<T>::~Allocator()
{
    for(unsigned int i = 0; i < prevTabs.size(); ++i)
    {
        delete [] prevTabs[i];
    }

    delete [] currentTab;
}

template <typename T>
void Allocator<T>::shiftTables()
{
    prevTabs.push_back(currentTab);
    currentTab      = new T[dim * tabSize];
    currentId       = 0;
}

template class Allocator<int>;

}  // namespace IsoSpec
