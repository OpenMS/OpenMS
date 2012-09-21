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
