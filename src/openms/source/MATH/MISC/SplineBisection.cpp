// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#include <OpenMS/MATH/MISC/SplineBisection.h>

#include <OpenMS/MATH/MISC/BSpline2d.h>
#include <OpenMS/MATH/MISC/CubicSpline2d.h>

namespace OpenMS::Math
{


    // explicit instantiation.
    template 
    void spline_bisection<BSpline2d>(const BSpline2d & peak_spline, 
        double const left_neighbor_mz,
        double const right_neighbor_mz,
        double & max_peak_mz,
        double & max_peak_int,
        double const threshold);


    // explicit instantiation.
    template 
    void spline_bisection<CubicSpline2d>(const CubicSpline2d & peak_spline, 
        double const left_neighbor_mz,
        double const right_neighbor_mz,
        double & max_peak_mz,
        double & max_peak_int,
        double const threshold);

} //OpenMS //Math
