// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License, or (at your option) any later version.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// --------------------------------------------------------------------------
// $Maintainer: Florian Zeller $
// $Authors: Lukas Mueller, Markus Mueller $
// --------------------------------------------------------------------------
//
///////////////////////////////////////////////////////////////////////////
//
//  written by Lukas N Mueller, 30.3.05
//  Lukas.Mueller@imsb.biol.ethz.ch
//  Group of Prof. Ruedi Aebersold, IMSB, ETH Hoenggerberg, Zurich
//
//  Ported to OpenMS by Florian Zeller, florian.zeller@bsse.ethz.ch
//  December 2010
//

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/LC_MS_XML_reader.h>

namespace OpenMS
{

// These values are overwritten by config
double LC_MS_XML_reader::TR_MIN = 0;
double LC_MS_XML_reader::TR_MAX = 0; // 180
double LC_MS_XML_reader::FEATURE_MZ_MIN = 0; // 200
double LC_MS_XML_reader::FEATURE_MZ_MAX = 0; //1800;
int LC_MS_XML_reader::FEATURE_CHRG_MIN = 0; //1;
int LC_MS_XML_reader::FEATURE_CHRG_MAX = 0;

}
