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


#include <iostream>
#include <stdlib.h>
#include "dirtyAllocator.h"


DirtyAllocator::DirtyAllocator(
    const int dim, const int tabSize
): tabSize(tabSize)
{
    cellSize        = sizeof(double) + sizeof(int) * dim;
    // Fix memory alignment problems for SPARC
    if(cellSize % sizeof(double) != 0)
    	cellSize += sizeof(double) - cellSize % sizeof(double);
    currentTab      = malloc( cellSize * tabSize );
    currentConf     = currentTab;
    endOfTablePtr = reinterpret_cast<char*>(currentTab) + cellSize*tabSize;
}


DirtyAllocator::~DirtyAllocator()
{
    for(unsigned int i = 0; i < prevTabs.size(); ++i) free(prevTabs[i]);
    free(currentTab);
}

void DirtyAllocator::shiftTables()
{
    prevTabs.push_back(currentTab);

    currentTab              = malloc( cellSize * tabSize );
    currentConf             = currentTab;
    endOfTablePtr   = reinterpret_cast<char*>(currentTab) + cellSize*tabSize;
}
