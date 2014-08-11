/************************************************************************
 * Copyright 2009 University Corporation for Atmospheric Research.
 * All rights reserved.
 *
 * Use of this code is subject to the standard BSD license:
 *
 *  http://www.opensource.org/licenses/bsd-license.html
 *
 * See the COPYRIGHT file in the source distribution for the license text,
 * or see this web page:
 *
 *  http://www.eol.ucar.edu/homes/granger/bspline/doc/
 *
 *************************************************************************/

// Instantiate the BSpline templates for type 

#include "BSplineBase.cpp"
#include "BSpline.cpp"

/// Instantiate BSplineBase for a library
template class  BSplineBase<double>;
template class  BSplineBase<float>;

/// Instantiate BSpline for a library
template class BSpline<double>;
template class BSpline<float>;
