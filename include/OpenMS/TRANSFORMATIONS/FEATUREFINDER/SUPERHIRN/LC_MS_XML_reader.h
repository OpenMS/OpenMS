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

// Flo: I keep this class because if its constants. No functionality though, just
// for maximum code compatiblity with the original superhirn.

#ifndef LC_MS_XML_READER_H
#define LC_MS_XML_READER_H

#include <OpenMS/CONCEPT/Types.h>

namespace OpenMS
{

class OPENMS_DLLAPI LC_MS_XML_reader{
  
  ////////////////////////////////////////////////
  // declaration of the public members:
  
public:
  
    static double TR_MIN;
  static double TR_MAX;
  static double FEATURE_MZ_MIN;
  static double FEATURE_MZ_MAX;
  static int FEATURE_CHRG_MIN;
  static int FEATURE_CHRG_MAX;
  // static double PEAK_SCORE_THERSHOLD;
//  static double PEAK_INTENSITY_THRESHOLD;
//  static bool EXTRACT_MONO_ISOTOPE_PROFILE;
//  static double SIGNAL_TO_NOISE_THERSHOLD;
//  static string DATA_STORAGE_XML_FORMAT_TYPE;
};

} // ns

#endif

    
