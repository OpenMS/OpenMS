// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg$
// $Authors: Timo Sachsenberg$
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

#include <OpenMS/ANALYSIS/NUXL/NuXLAnnotateAndLocate.h>
#include <OpenMS/ANALYSIS/NUXL/NuXLAnnotatedHit.h>
#include <OpenMS/ANALYSIS/NUXL/NuXLParameterParsing.h>
#include <OpenMS/ANALYSIS/NUXL/NuXLPresets.h>
#include <OpenMS/ANALYSIS/NUXL/NuXLModificationsGenerator.h>
#include <OpenMS/ANALYSIS/NUXL/NuXLPreprocessSpectra.h>
#include <OpenMS/ANALYSIS/ID/PrecursorPurity.h>
#include <OpenMS/CHEMISTRY/ModifiedPeptideGenerator.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/Peak1D.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/DATASTRUCTURES/StringView.h>

using namespace std;
using namespace OpenMS;

START_TEST(NuXLAnnotateAndLocate, "$Id$")

/////////////////////////////////////////////////////////////

// Note: NuXLAnnotateAndLocate only contains static methods, so no constructor/destructor tests needed

START_SECTION((static void annotateAndLocate_(const PeakMap& exp, std::vector<std::vector<NuXLAnnotatedHit>>& annotated_hits, const NuXLModificationMassesResult& mm, const ModifiedPeptideGenerator::MapToResidueType& fixed_modifications, const ModifiedPeptideGenerator::MapToResidueType& variable_modifications, Size max_variable_mods_per_peptide, double fragment_mass_tolerance, bool fragment_mass_tolerance_unit_ppm, const NuXLParameterParsing::PrecursorsToMS2Adducts& all_feasible_adducts)))
{
  // Create test data for the annotateAndLocate_ method
  PeakMap exp;
  
  // Create MS2 spectrum
  MSSpectrum spectrum;
  
  // Set spectrum metadata
  spectrum.setMSLevel(2);
  spectrum.setRT(26.236534118652301); // RT in seconds
  spectrum.setNativeID("scan=0");
  
  // Set precursor information
  Precursor precursor;
  precursor.setMZ(678.284423828125);
  precursor.setCharge(2);
  precursor.setIntensity(0.0); // Not specified in MGF
  
  vector<Precursor> precursors;
  precursors.push_back(precursor);
  spectrum.setPrecursors(precursors);
  
  // Define peak data (m/z, intensity pairs)
  vector<pair<double, double>> peak_data = {
      {110.071525573730469, 6.253766e04},
      {112.050888061523438, 3.577897e04},
      {112.086700439453125, 1.105394e04},
      {114.056053161621094, 2988.362305},
      {115.086708068847656, 3076.342773},
      {116.070907592773438, 3446.308594},
      {129.102294921875, 1.756605e04},
      {130.032073974609375, 4686.557129},
      {130.08587646484375, 5617.307129},
      {136.0618896484375, 2.164918e05},
      {136.074996948242188, 2.029108e04},
      {137.047286987304688, 3687.551514},
      {137.065200805664063, 6387.198242},
      {147.1129150390625, 1.92404e04},
      {148.042221069335938, 3363.700928},
      {152.056503295898438, 1.203792e04},
      {158.092330932617188, 1.480457e04},
      {161.092391967773438, 3599.997803},
      {164.99444580078125, 3784.263916},
      {169.13360595703125, 4.185202e04},
      {173.092025756835938, 1.233966e04},
      {175.071029663085938, 1.580166e04},
      {175.118759155273438, 8.477893e04},
      {179.045318603515625, 1.823629e04},
      {189.0687255859375, 1.301581e04},
      {192.101638793945313, 1.350166e04},
      {197.128097534179688, 1.393258e04},
      {198.168533325195313, 1.331057e04},
      {198.191116333007813, 2.554027e04},
      {204.134170532226563, 2.405196e04},
      {205.065078735351563, 3374.075928},
      {213.014556884765625, 4348.111328},
      {218.151031494140625, 3751.484131},
      {227.066131591796875, 1.017478e05},
      {228.071426391601563, 1.20593e04},
      {229.127822875976563, 3365.893799},
      {232.13983154296875, 1.504747e04},
      {233.105987548828125, 5794.908691},
      {233.131546020507813, 1.410256e05},
      {237.123031616210938, 1.239832e04},
      {239.113800048828125, 3990.227051},
      {240.134185791015625, 5963.953613},
      {246.155319213867188, 1.500475e04},
      {250.093399047851563, 5158.16748},
      {251.101776123046875, 1.641352e04},
      {261.126708984375, 8.40667e04},
      {268.164642333984375, 1.894302e04},
      {284.12322998046875, 3232.119629},
      {286.150421142578125, 1.869846e04},
      {288.20086669921875, 3759.060547},
      {303.177154541015625, 1.501742e05},
      {304.178009033203125, 1.185856e04},
      {306.1414794921875, 3565.862061},
      {312.047698974609375, 5508.858398},
      {317.217926025390625, 1.34463e04},
      {330.05950927734375, 4.907482e04},
      {332.163787841796875, 5.271074e04},
      {345.22412109375, 2.51371e04},
      {348.0694580078125, 3864.993897},
      {359.13421630859375, 4295.918457},
      {362.1336669921875, 1.616536e04},
      {363.131317138671875, 1.351985e04},
      {368.21673583984375, 5.284671e04},
      {368.71990966796875, 2.359747e04},
      {377.221343994140625, 7.732263e04},
      {377.719268798828125, 2.048313e04},
      {399.233306884765625, 1.637699e04},
      {415.200531005859375, 5554.976074},
      {416.16107177734375, 4973.861328},
      {416.2611083984375, 1.733192e05},
      {417.262054443359375, 2.281113e04},
      {418.123992919921875, 1.756198e04},
      {423.222259521484375, 4217.286133},
      {424.219085693359375, 1.986142e04},
      {435.182647705078125, 5197.253906},
      {459.148712158203125, 4979.168457},
      {460.137542724609375, 3691.867432},
      {466.247894287109375, 3510.215576},
      {477.163330078125, 2.09019e04},
      {484.256866455078125, 2.417197e04},
      {488.129302978515625, 4284.745117},
      {499.19189453125, 4154.002441},
      {505.15753173828125, 2.958442e04},
      {507.261810302734375, 4367.555176},
      {510.177032470703125, 4766.870606},
      {514.7666015625, 4061.070557},
      {515.2484130859375, 4104.623047},
      {515.33038330078125, 1.04829e05},
      {516.34088134765625, 1.926769e04},
      {517.205078125, 5048.685547},
      {519.75439453125, 5256.198242},
      {527.185791015625, 1.428935e04},
      {528.1842041015625, 1.262405e04},
      {541.75848388671875, 1.145473e04},
      {544.24188232421875, 1.113731e04},
      {550.15606689453125, 3788.122559},
      {554.76214599609375, 5055.256836},
      {555.177978515625, 5592.391113},
      {562.26092529296875, 1.931875e04},
      {562.775390625, 1.217805e04},
      {563.275390625, 3.180423e04},
      {563.7735595703125, 1.516493e04},
      {567.764404296875, 6875.622559},
      {568.26849365234375, 6272.043457},
      {571.7821044921875, 5.705805e04},
      {572.27984619140625, 6.670002e04},
      {572.78509521484375, 3.137877e04},
      {576.77349853515625, 3.774367e04},
      {577.27886962890625, 2.628601e04},
      {580.26910400390625, 1.242485e04},
      {583.31561279296875, 4133.029297},
      {585.777587890625, 4305.655762},
      {589.7821044921875, 1.17204e04},
      {590.28338623046875, 4917.491699},
      {598.28643798828125, 2.424438e04},
      {598.37109375, 1.423503e04},
      {602.77374267578125, 1.171851e04},
      {603.2696533203125, 5092.021484},
      {611.28338623046875, 1.388459e04},
      {611.7808837890625, 3.834531e04},
      {612.27587890625, 1.353678e04},
      {614.22186279296875, 5917.741699},
      {616.3775634765625, 2.142273e05},
      {617.3787841796875, 4.907813e04},
      {620.29071044921875, 7.968598e04},
      {620.7886962890625, 7.448647e04},
      {621.28009033203125, 1.615241e04},
      {624.2008056640625, 2.882028e04},
      {625.20501708984375, 6245.141602},
      {626.36175537109375, 5554.005859},
      {629.2969970703125, 3.109704e05},
      {629.7974853515625, 1.909436e05},
      {638.30474853515625, 1.256815e04},
      {642.2138671875, 2.738232e04},
      {643.2093505859375, 5399.541504},
      {648.33551025390625, 1.2203e04},
      {656.23358154296875, 5101.35791},
      {660.23907470703125, 1.513047e04},
      {660.77899169921875, 4942.533203},
      {669.27935791015625, 1.950834e04},
      {669.778564453125, 1.390851e04},
      {678.237548828125, 4.527829e04},
      {696.25799560546875, 5.060615e04},
      {722.1842041015625, 6325.642578},
      {725.2520751953125, 2.302715e04},
      {735.421142578125, 1.520849e04},
      {740.19207763671875, 5457.153809},
      {743.26336669921875, 2.377473e04},
      {753.4393310546875, 1.013214e05},
      {754.43768310546875, 3.028721e04},
      {777.3187255859375, 4896.175293},
      {795.32861328125, 1.31911e04},
      {798.3021240234375, 5863.57666},
      {813.4117431640625, 4935.314941},
      {814.327880859375, 6533.92627},
      {824.3211669921875, 1.895101e04},
      {825.31719970703125, 1.169664e04},
      {834.40716552734375, 4818.745117},
      {842.3365478515625, 2.336112e04},
      {843.3336181640625, 5359.795898},
      {884.44525146484375, 6872.038086},
      {898.3970947265625, 1.272029e04},
      {909.42108154296875, 5247.175293},
      {916.491455078125, 5883.089844},
      {921.43280029296875, 4168.393555},
      {927.43896484375, 1.906384e04},
      {937.42437744140625, 5556.257813},
      {955.418212890625, 2.53536e04},
      {980.47265625, 4401.492188},
      {985.42645263671875, 1.125866e04},
      {998.47479248046875, 1.84024e04},
      {1025.4136962890625, 1.467806e04},
      {1031.5264892578125, 6.322954e04},
      {1032.5228271484375, 2.714641e04},
      {1078.453857421875, 5037.069824},
      {1083.4879150390625, 1.257548e04},
      {1096.454833984375, 1.260347e04},
      {1143.546875, 1.702259e04},
      {1144.5631103515625, 6057.237793}
  };
  
  // Add peaks to spectrum
  spectrum.reserve(peak_data.size());
  for (const auto& peak : peak_data) {
      spectrum.push_back(Peak1D(peak.first, peak.second));
  }
  
  // Sort peaks by m/z (recommended for OpenMS spectra)
  spectrum.sortByPosition();
  
  exp.addSpectrum(spectrum);
  
  // Preprocess spectra analogous to OpenNuXL using the newly added method
  // This will properly annotate charge arrays and perform other preprocessing steps
  std::map<String, PrecursorPurity::PurityScores> purities; // Empty purities map for test
  double window_size = 100.0; // Default window size from OpenNuXL
  Size peak_count = 6; // Default peak count from OpenNuXL
  
  NuXLPreprocessSpectra::preprocessSpectra(exp,
                                            false, // no single charge conversion
                                            true,  // annotate charge
                                            window_size,
                                            peak_count,
                                            purities);
  
  // Create empty annotated hits vector
  vector<vector<NuXLAnnotatedHit>> annotated_hits;
  annotated_hits.resize(1); // One spectrum
  
  // Retrieve preset for RNA-UV (U) to get default PrecursorToMS2Adducts
  StringList modifications;
  StringList fragment_adducts;
  String can_cross_link;
  StringList target_nucleotides;
  StringList mappings;

  NuXLPresets::getPresets("RNA-UV (U)", target_nucleotides, mappings, modifications, fragment_adducts, can_cross_link);

  // Convert string to set
  set<char> can_xl;
  for (const auto& c : can_cross_link) { can_xl.insert(c); }

  // Create NuXLModificationMassesResult
  NuXLModificationMassesResult mm = NuXLModificationsGenerator::initModificationMassesNA(
            target_nucleotides,
            StringList(),
            can_xl,
            mappings,
            modifications,
            "",
            false,
            2);

  mm.formula2mass[""] = 0; // insert "null" modification otherwise peptides without NA will not be searched
  mm.mod_combinations[""].insert("none");
  
  // Create empty modification maps
  ModifiedPeptideGenerator::MapToResidueType fixed_modifications;
  ModifiedPeptideGenerator::MapToResidueType variable_modifications;
  
  // Get nucleotide to fragment adducts mapping
  NuXLParameterParsing::NucleotideToFragmentAdductMap nucleotide_to_fragment_adducts =
    NuXLParameterParsing::getTargetNucleotideToFragmentAdducts(fragment_adducts);
  
  // Get all feasible fragment adducts from all possible precursor adducts
  NuXLParameterParsing::PrecursorsToMS2Adducts all_feasible_adducts =
    NuXLParameterParsing::getAllFeasibleFragmentAdducts(mm, nucleotide_to_fragment_adducts, can_xl, true, true);
  
  // Test parameters
  Size max_variable_mods_per_peptide = 2;
  double fragment_mass_tolerance = 20.0; // 20 ppm
  bool fragment_mass_tolerance_unit_ppm = true;
  
  // Create an annotated hit with peptide sequence DYHTVLGAR and precursor adduct U-H2O1
  NuXLAnnotatedHit hit;
  const String seq = "DYHTVLGAR";
  hit.sequence = StringView(seq);
  hit.peptide_mod_index = 0; // No peptide modifications
  hit.cross_linked_nucleotide = 'U';
  
  // Find the index for "U-H2O1" in the modification combinations
  // "U-H2O1" corresponds to formula "C9H11N2O8P1"
  Size na_mod_index = 0;
  for (const auto& pair : mm.mod_combinations) {
    if (pair.first == "C9H11N2O8P1" && pair.second.find("U-H2O1") != pair.second.end()) {
      break;
    }
    na_mod_index++;
  }
  hit.NA_mod_index = na_mod_index;
  
  // Set some basic scores for a realistic hit
  hit.score = 10.0;
  hit.mass_error_p = 0.5;
  hit.total_loss_score = 5.0;
  hit.total_MIC = 0.3;
  
  // Add the hit to the annotated hits
  annotated_hits[0].push_back(hit);
  
  // Test that the method can handle a real annotated hit after preprocessing
  NuXLAnnotateAndLocate::annotateAndLocate_(
    exp,
    annotated_hits,
    mm,
    fixed_modifications,
    variable_modifications,
    max_variable_mods_per_peptide,
    fragment_mass_tolerance,
    fragment_mass_tolerance_unit_ppm,
    all_feasible_adducts
  );

    // Test localization and annotation results (dummy expected values for now)
  TEST_EQUAL(annotated_hits[0].size(), 1)
  const NuXLAnnotatedHit& processed_hit = annotated_hits[0][0];

    // Expected best_localization (dummy value - to be replaced with correct one)
  String expected_best_localization = "DyHTVLGAR"; // or specific localization pattern
  TEST_EQUAL(processed_hit.best_localization, expected_best_localization)

  String expected_localization_scores = "0,73.58,0,0,0,0,0,0,0";
  TEST_EQUAL(processed_hit.localization_scores, expected_localization_scores)
  
  // Expected best_localization_score
  double expected_best_localization_score = 0.73578;
  TEST_REAL_SIMILAR(processed_hit.best_localization_score, expected_best_localization_score)
  
  // Expected best_localization_position
  int expected_best_localization_position = 1;
  TEST_EQUAL(processed_hit.best_localization_position, expected_best_localization_position)
  
  // Test that fragment annotations were generated with expected values
  TEST_EQUAL(processed_hit.fragment_annotations.size(), 30)
  
  // Test all fragment annotations with their expected values
  TEST_EQUAL(processed_hit.fragment_annotations[0].annotation, "y1+")
  TEST_REAL_SIMILAR(processed_hit.fragment_annotations[0].mz, 175.118759155273)
  TEST_REAL_SIMILAR(processed_hit.fragment_annotations[0].intensity, 0.168911263346672)
  
  TEST_EQUAL(processed_hit.fragment_annotations[1].annotation, "y3+")
  TEST_REAL_SIMILAR(processed_hit.fragment_annotations[1].mz, 303.177154541016)
  TEST_REAL_SIMILAR(processed_hit.fragment_annotations[1].intensity, 0.322829753160477)
  
  TEST_EQUAL(processed_hit.fragment_annotations[2].annotation, "y4+")
  TEST_REAL_SIMILAR(processed_hit.fragment_annotations[2].mz, 416.261108398438)
  TEST_REAL_SIMILAR(processed_hit.fragment_annotations[2].intensity, 0.390764802694321)
  
  TEST_EQUAL(processed_hit.fragment_annotations[3].annotation, "y5+")
  TEST_REAL_SIMILAR(processed_hit.fragment_annotations[3].mz, 515.330383300781)
  TEST_REAL_SIMILAR(processed_hit.fragment_annotations[3].intensity, 0.247246921062469)
  
  TEST_EQUAL(processed_hit.fragment_annotations[4].annotation, "y6+")
  TEST_REAL_SIMILAR(processed_hit.fragment_annotations[4].mz, 616.377563476562)
  TEST_REAL_SIMILAR(processed_hit.fragment_annotations[4].intensity, 0.524602711200714)
  
  TEST_EQUAL(processed_hit.fragment_annotations[5].annotation, "y7++")
  TEST_REAL_SIMILAR(processed_hit.fragment_annotations[5].mz, 377.221343994141)
  TEST_REAL_SIMILAR(processed_hit.fragment_annotations[5].intensity, 0.194865584373474)
  
  TEST_EQUAL(processed_hit.fragment_annotations[6].annotation, "y7+")
  TEST_REAL_SIMILAR(processed_hit.fragment_annotations[6].mz, 753.439331054688)
  TEST_REAL_SIMILAR(processed_hit.fragment_annotations[6].intensity, 0.262213468551636)
  
  TEST_EQUAL(processed_hit.fragment_annotations[7].annotation, "y8+")
  TEST_REAL_SIMILAR(processed_hit.fragment_annotations[7].mz, 916.491455078125)
  TEST_REAL_SIMILAR(processed_hit.fragment_annotations[7].intensity, 0.0117213102057576)
  
  TEST_EQUAL(processed_hit.fragment_annotations[8].annotation, "iH+")
  TEST_REAL_SIMILAR(processed_hit.fragment_annotations[8].mz, 110.07152557373)
  TEST_REAL_SIMILAR(processed_hit.fragment_annotations[8].intensity, 0.124598354101181)
  
  TEST_EQUAL(processed_hit.fragment_annotations[9].annotation, "y7-H2O1++")
  TEST_REAL_SIMILAR(processed_hit.fragment_annotations[9].mz, 368.216735839844)
  TEST_REAL_SIMILAR(processed_hit.fragment_annotations[9].intensity, 0.152305334806442)
  
  TEST_EQUAL(processed_hit.fragment_annotations[10].annotation, "y7-H2O1+")
  TEST_REAL_SIMILAR(processed_hit.fragment_annotations[10].mz, 735.421142578125)
  TEST_REAL_SIMILAR(processed_hit.fragment_annotations[10].intensity, 0.0303009878844023)
  
  TEST_EQUAL(processed_hit.fragment_annotations[11].annotation, "b2+U-H3PO4+")
  TEST_REAL_SIMILAR(processed_hit.fragment_annotations[11].mz, 505.157531738281)
  TEST_REAL_SIMILAR(processed_hit.fragment_annotations[11].intensity, 0.0589432045817375)
  
  TEST_EQUAL(processed_hit.fragment_annotations[12].annotation, "b4+U-H3PO4+")
  TEST_REAL_SIMILAR(processed_hit.fragment_annotations[12].mz, 743.263366699219)
  TEST_REAL_SIMILAR(processed_hit.fragment_annotations[12].intensity, 0.0473681353032589)
  
  TEST_EQUAL(processed_hit.fragment_annotations[13].annotation, "b5+U-H3PO4+")
  TEST_REAL_SIMILAR(processed_hit.fragment_annotations[13].mz, 842.336547851562)
  TEST_REAL_SIMILAR(processed_hit.fragment_annotations[13].intensity, 0.0572227798402309)
  
  TEST_EQUAL(processed_hit.fragment_annotations[14].annotation, "b6+U-H3PO4+")
  TEST_REAL_SIMILAR(processed_hit.fragment_annotations[14].mz, 955.418212890625)
  TEST_REAL_SIMILAR(processed_hit.fragment_annotations[14].intensity, 0.0505138337612152)
  
  TEST_EQUAL(processed_hit.fragment_annotations[15].annotation, "b7+U'+")
  TEST_REAL_SIMILAR(processed_hit.fragment_annotations[15].mz, 898.397094726562)
  TEST_REAL_SIMILAR(processed_hit.fragment_annotations[15].intensity, 0.0253435652703047)
  
  TEST_EQUAL(processed_hit.fragment_annotations[16].annotation, "b8+C3O+")
  TEST_REAL_SIMILAR(processed_hit.fragment_annotations[16].mz, 909.421081542969)
  TEST_REAL_SIMILAR(processed_hit.fragment_annotations[16].intensity, 0.0104543315246701)
  
  TEST_EQUAL(processed_hit.fragment_annotations[17].annotation, "b8+U-H3PO4+")
  TEST_REAL_SIMILAR(processed_hit.fragment_annotations[17].mz, 1083.48791503906)
  TEST_REAL_SIMILAR(processed_hit.fragment_annotations[17].intensity, 0.0250550508499146)
  
  TEST_EQUAL(processed_hit.fragment_annotations[18].annotation, "y8+U-H3PO4++")
  TEST_REAL_SIMILAR(processed_hit.fragment_annotations[18].mz, 571.782104492188)
  TEST_REAL_SIMILAR(processed_hit.fragment_annotations[18].intensity, 0.309090495109558)
  
  TEST_EQUAL(processed_hit.fragment_annotations[19].annotation, "a2+U-H3PO4+")
  TEST_REAL_SIMILAR(processed_hit.fragment_annotations[19].mz, 477.163330078125)
  TEST_REAL_SIMILAR(processed_hit.fragment_annotations[19].intensity, 0.041644386947155)
  
  TEST_EQUAL(processed_hit.fragment_annotations[20].annotation, "a5+U-H3PO4+")
  TEST_REAL_SIMILAR(processed_hit.fragment_annotations[20].mz, 814.327880859375)
  TEST_REAL_SIMILAR(processed_hit.fragment_annotations[20].intensity, 0.0130180194973946)
  
  TEST_EQUAL(processed_hit.fragment_annotations[21].annotation, "a6+U-H3PO4+")
  TEST_REAL_SIMILAR(processed_hit.fragment_annotations[21].mz, 927.43896484375)
  TEST_REAL_SIMILAR(processed_hit.fragment_annotations[21].intensity, 0.0379822850227356)
  
  TEST_EQUAL(processed_hit.fragment_annotations[22].annotation, "iY+U-H3PO4++")
  TEST_REAL_SIMILAR(processed_hit.fragment_annotations[22].mz, 362.133666992188)
  TEST_REAL_SIMILAR(processed_hit.fragment_annotations[22].intensity, 0.0591440200805664)
  
  TEST_EQUAL(processed_hit.fragment_annotations[23].annotation, "MI:C'+")
  TEST_REAL_SIMILAR(processed_hit.fragment_annotations[23].mz, 112.050888061523)
  TEST_REAL_SIMILAR(processed_hit.fragment_annotations[23].intensity, 0.0712850615382195)
  
  TEST_EQUAL(processed_hit.fragment_annotations[24].annotation, "MI:A'+")
  TEST_REAL_SIMILAR(processed_hit.fragment_annotations[24].mz, 136.061889648438)
  TEST_REAL_SIMILAR(processed_hit.fragment_annotations[24].intensity, 0.444058150053024)
  
  TEST_EQUAL(processed_hit.fragment_annotations[25].annotation, "MI:U-H3PO4+")
  TEST_REAL_SIMILAR(processed_hit.fragment_annotations[25].mz, 227.066131591797)
  TEST_REAL_SIMILAR(processed_hit.fragment_annotations[25].intensity, 0.226746201515198)
  
  TEST_EQUAL(processed_hit.fragment_annotations[26].annotation, "[M+H]+")
  TEST_REAL_SIMILAR(processed_hit.fragment_annotations[26].mz, 1031.52648925781)
  TEST_REAL_SIMILAR(processed_hit.fragment_annotations[26].intensity, 0.180062621831894)
  
  TEST_EQUAL(processed_hit.fragment_annotations[27].annotation, "[M+2H-H2O]+U-H3PO4++")
  TEST_REAL_SIMILAR(processed_hit.fragment_annotations[27].mz, 620.290710449219)
  TEST_REAL_SIMILAR(processed_hit.fragment_annotations[27].intensity, 0.307169020175934)
  
  TEST_EQUAL(processed_hit.fragment_annotations[28].annotation, "[M+2H]+U-H3PO4++")
  TEST_REAL_SIMILAR(processed_hit.fragment_annotations[28].mz, 629.296997070312)
  TEST_REAL_SIMILAR(processed_hit.fragment_annotations[28].intensity, 1)
  
  TEST_EQUAL(processed_hit.fragment_annotations[29].annotation, "[M+2H+U-H2O1]++")
  TEST_REAL_SIMILAR(processed_hit.fragment_annotations[29].mz, 669.279357910156)
  TEST_REAL_SIMILAR(processed_hit.fragment_annotations[29].intensity, 0.0665788426995277)
  
}
END_SECTION

END_TEST