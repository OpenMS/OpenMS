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

#include "platform.h"

#if !ISOSPEC_BUILDING_R

#if !ISOSPEC_GOT_SYSTEM_MMAN && ISOSPEC_GOT_MMAN
    #include "mman.cpp"  // NOLINT(build/include)
#endif

// A poor-man's replacement for LTO. We're small enough that we can do that. And
// ignore cpplint's complaints about it.

#include "allocator.cpp"        // NOLINT(build/include)
#include "dirtyAllocator.cpp"   // NOLINT(build/include)
#include "isoSpec++.cpp"        // NOLINT(build/include)
#include "isoMath.cpp"          // NOLINT(build/include)
#include "marginalTrek++.cpp"   // NOLINT(build/include)
#include "operators.cpp"        // NOLINT(build/include)
#include "element_tables.cpp"   // NOLINT(build/include)
#include "fasta.cpp"            // NOLINT(build/include)
#include "cwrapper.cpp"         // NOLINT(build/include)
#include "fixedEnvelopes.cpp"   // NOLINT(build/include)
#include "misc.cpp"             // NOLINT(build/include)

#endif
