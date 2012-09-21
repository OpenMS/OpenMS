///////////////////////////////////////////////////////////////////////////
//
//  written by Lukas N Mueller, 30.3.05
//  Lukas.Mueller@imsb.biol.ethz.ch
//  Group of Prof. Ruedi Aebersold, IMSB, ETH Hoenggerberg, Zurich
//
//  Ported to OpenMS by Florian Zeller, florian.zeller@bsse.ethz.ch
//  December 2010
//
//
// **********************************************************************//
// CLASS LC_MS_XML_reader:
//
// variable description:
//
//
// function description:
//
//
// **********************************************************************//



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

    
