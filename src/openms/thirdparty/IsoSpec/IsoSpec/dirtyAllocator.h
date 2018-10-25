/*
 *   Copyright (C) 2015-2016 Mateusz Łącki and Michał Startek.
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


#ifndef DIRTY_ALLOCATOR_HPP
#define DIRTY_ALLOCATOR_HPP

#include <vector>
#include <iostream>
#include <string.h>

class DirtyAllocator{
private:
    void*   currentTab;
    void*   currentConf;
    void*   endOfTablePtr;
    const int       tabSize;
    int     cellSize;
    std::vector<void*>  prevTabs;
public:
    DirtyAllocator(const int dim, const int tabSize = 10000);
    ~DirtyAllocator();

    void shiftTables();

    inline void* newConf()
    {
        if (currentConf >= endOfTablePtr)
        {
            shiftTables();
        }

        void*  ret = currentConf;
        currentConf = reinterpret_cast<char*>(currentConf) + cellSize;

        return ret;
    }

    inline void* makeCopy(const void* conf)
    {
        void* currentPlace = newConf();

        memcpy(currentPlace, conf, cellSize);

        return currentPlace;
    }

    inline void* makeExternalCopy(const void* conf)
    {
        void* res = malloc(cellSize);

        memcpy(res, conf, cellSize);

        return res;
    }
};

#endif
