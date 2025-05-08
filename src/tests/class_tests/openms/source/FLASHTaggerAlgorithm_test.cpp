// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Jaekwan Kim$
// $Authors: Jaekwan Kim$
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHTaggerAlgorithm.h>
#include <OpenMS/FORMAT/MzMLFile.h>
//#include <OpenMS/FORMAT/MzMLFile.h>
///////////////////////////


using namespace OpenMS;
using namespace std;

START_TEST(FLASHTaggerAlgorithm, "$Id$")

// /////////////////////////////////////////////////////////////
// /////////////////////////////////////////////////////////////

// FLASHTaggerAlgorithm* ptr = 0;
// FLASHTaggerAlgorithm* null_ptr = 0;
// START_SECTION(FLASHTaggerAlgorithm())
// {
//   ptr = new FLASHTaggerAlgorithm();
//   TEST_NOT_EQUAL(ptr, null_ptr)
// }
// END_SECTION

// START_SECTION(~FLASHTaggerAlgorithm())
// {
//   delete ptr;
// }
// END_SECTION

// PeakMap map;
// MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("SIGMAmAb_ECDHCD_nativeTD_1_deconv.mzML"), map);

// DeconvolvedSpectrum deconvolved_spectrum;
// double tol(10);
// // run FLASHDeconvAlgorithm here!


// // collect statistics for information
// for (int index = 0; index < map.size(); index++)
// {

//   auto spec = map[index];
//   if (spec.getMSLevel() != 2) continue;
//   int scan = FLASHDeconvAlgorithm::getScanNumber(map, index);
//   DeconvolvedSpectrum dspec(scan);
//   dspec.setOriginalSpectrum(spec);
//   String deconv_meta_str = spec.getMetaValue("DeconvMassInfo").toString();
//   int tol_loc_s = deconv_meta_str.find("tol=") + 4;
//   int tol_loc_e = deconv_meta_str.find(";", tol_loc_s);
//   tol = stod(deconv_meta_str.substr(tol_loc_s, tol_loc_e - tol_loc_s));

//   int q_loc_s = deconv_meta_str.find("qscore=") + 7;
//   int q_loc_e = deconv_meta_str.find(";", q_loc_s);
//   auto q_str = deconv_meta_str.substr(q_loc_s, q_loc_e - q_loc_s);
//   Size pos = 0;
//   std::vector<double> qscores;
//   while (true)
//   {
//     Size pos_t = q_str.find(",", pos);
//     if (pos_t == String::npos) break;
//     auto token = q_str.substr(pos, pos_t - pos);
//     qscores.push_back(stod(token));
//     pos = pos_t + 1;
//   }

//   int s_loc_s = deconv_meta_str.find("snr=") + 4;
//   int s_loc_e = deconv_meta_str.find(";", s_loc_s);
//   auto s_str = deconv_meta_str.substr(s_loc_s, s_loc_e - s_loc_s);
//   pos = 0;
//   std::vector<float> snrs;
//   while (true)
//   {
//     Size pos_t = s_str.find(",", pos);
//     if (pos_t == String::npos) break;
//     auto token = s_str.substr(pos, pos_t - pos);
//     snrs.push_back(stof(token));
//     pos = pos_t + 1;
//   }

//   for (int i = 0; i < spec.size(); i++)
//   {
//     PeakGroup peak;
//     peak.setQscore(qscores[i]);
//     peak.setSNR(snrs[i]);
//     peak.setMonoisotopicMass(spec[i].getMZ());
//     peak.setScanNumber(scan);
//     dspec.push_back(peak);
//   }
//   dspec.sort();
//   deconvolved_spectrum = dspec;
// }
// // Run tagger
// FLASHTaggerAlgorithm tagger;




// START_SECTION(run())
// {
//   FASTAFile fasta_file;
//   std::vector<FASTAFile::FASTAEntry> fasta_entry;
//   fasta_file.load("C:\\Users\\qlcsk\\Desktop\\jkvision\\openms\\pyFLASHDeconv\\uniprot-mg1655-filtered-reviewed_yes.fasta", fasta_entry);
//   tagger.run(deconvolved_spectrum, tol, fasta_entry);
//   std::vector<FLASHHelperClasses::Tag> tag;
//   tagger.getTags(true, tag);
//   TEST_EQUAL(tag.size(),114)
// }
// END_SECTION







// /////////////////////////////////////////////////////////////
// /////////////////////////////////////////////////////////////
END_TEST