// David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2022
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt
// https://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// Version: 6.0.2022.01.06

#pragma once

#include <Mathematics/Math.h>

namespace gte
{
    // Test-intersection queries.

    template <typename Real, typename Type0, typename Type1>
    class TIQuery
    {
    public:
        struct Result
        {
            // A TIQuery-base class B must define a B::Result struct with
            // member 'bool intersect'.  A TIQuery-derived class D must also
            // derive a D::Result from B:Result but may have no members.  The
            // member 'intersect' is 'true' iff the primitives intersect.  The
            // operator() is const for conceptual constness, but derived
            // classes can use internal data to support the queries and tag
            // that data with the mutable modifier.
        };

        Result operator()(Type0 const& primitive0, Type1 const& primitive1) const;
    };
}
