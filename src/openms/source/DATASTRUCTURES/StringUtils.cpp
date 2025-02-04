// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg, Chris Bielow $
// $Authors: Marc Sturm, Stephan Aiche, Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/DATASTRUCTURES/StringUtils.h>

namespace OpenMS
{

    boost::spirit::qi::real_parser<double, StringUtilsHelper::real_policies_NANfixed_<double> > StringUtilsHelper::parse_double_ = boost::spirit::qi::real_parser<double, real_policies_NANfixed_<double> >();
    boost::spirit::qi::real_parser<float, StringUtilsHelper::real_policies_NANfixed_<float> > StringUtilsHelper::parse_float_ = boost::spirit::qi::real_parser<float, real_policies_NANfixed_<float> >();
 
}
