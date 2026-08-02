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

#include <vector>
#include <cstring>

namespace IsoSpec
{

class DirtyAllocator
{
 private:
    void*   currentTab;
    void*   currentConf;
    void*   endOfTablePtr;
    const int       tabSize;
    int     cellSize;
    std::vector<void*>  prevTabs;

 public:
    explicit DirtyAllocator(const int dim, const int tabSize = 10000);
    ~DirtyAllocator();

    DirtyAllocator(const DirtyAllocator& other) = delete;
    DirtyAllocator& operator=(const DirtyAllocator& other) = delete;

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
};

}  // namespace IsoSpec
