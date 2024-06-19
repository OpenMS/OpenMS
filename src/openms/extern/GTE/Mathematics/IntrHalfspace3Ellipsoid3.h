// David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2022
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt
// https://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// Version: 6.0.2022.01.06

#pragma once

#include <Mathematics/TIQuery.h>
#include <Mathematics/Halfspace.h>
#include <Mathematics/Hyperellipsoid.h>
#include <Mathematics/Matrix3x3.h>

// Queries for intersection of objects with halfspaces.  These are useful for
// containment testing, object culling, and clipping.

namespace gte
{
    template <typename T>
    class TIQuery<T, Halfspace3<T>, Ellipsoid3<T>>
    {
    public:
        struct Result
        {
            Result()
                :
                intersect(false)
            {
            }

            bool intersect;
        };

        Result operator()(Halfspace3<T> const& halfspace, Ellipsoid3<T> const& ellipsoid)
        {
            // Project the ellipsoid onto the normal line.  The plane of the
            // halfspace occurs at the origin (zero) of the normal line.
            Result result{};
            Matrix3x3<T> MInverse;
            ellipsoid.GetMInverse(MInverse);
            T discr = Dot(halfspace.normal, MInverse * halfspace.normal);
            T extent = std::sqrt(std::max(discr, (T)0));
            T center = Dot(halfspace.normal, ellipsoid.center) - halfspace.constant;
            T tmax = center + extent;

            // The ellipsoid and halfspace intersect when the projection
            // interval maximum is nonnegative.
            result.intersect = (tmax >= (T)0);
            return result;
        }
    };
}
