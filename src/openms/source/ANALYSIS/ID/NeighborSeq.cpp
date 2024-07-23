// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Sven Nahnsen $
// $Authors: Sven Nahnsen, Andreas Bertsch, Chris Bielow, Philipp Wang $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/ANALYSIS/ID/NeighborSeq.h>
#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGenerator.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <cmath>



using namespace OpenMS;
using namespace std;




// Function to generate the theoretical spectrum for a given peptide sequence
MSSpectrum NeighborSeq::generateSpectrum(const String& peptide_sequence)
{
AASequence peptide = AASequence::fromString(peptide_sequence);
TheoreticalSpectrumGenerator spec_gen;
// Set parameters for the spectrum generator
Param params;
params.setValue("add_b_ions", "true");
params.setValue("add_y_ions", "true");
spec_gen.setParameters(params);
// Generate the spectrum and store it in an MSSpectrum object
MSSpectrum spectrum;
spec_gen.getSpectrum(spectrum, peptide, 1, 1);
return spectrum;
}

// Function to compare two spectra and determine if they are similar
bool NeighborSeq::compareSpectra(const MSSpectrum& spec1, const MSSpectrum& spec2, const double& min_shared_ion_fraction, const double& mz_bin_size)
{
  std::map<int, int> bins1, bins2;

  // Lambda function to bin a spectrum into a map
  auto bin_spectrum = [&](const MSSpectrum& spec, std::map<int, int>& bins) {
    for (const auto& peak : spec)
      bins[static_cast<int>(peak.getMZ() / mz_bin_size)]++;
  };

  // Bin both spectra
  bin_spectrum(spec1, bins1);
  bin_spectrum(spec2, bins2);

  // Extract bin keys as vectors
  std::vector<int> vec_bins1, vec_bins2;
  for (const auto& bin : bins1)
    vec_bins1.push_back(bin.first);
  for (const auto& bin : bins2)
    vec_bins2.push_back(bin.first);

  // Calculate the intersection of the bin vectors
  std::vector<int> intersection;
  std::set_intersection(vec_bins1.begin(), vec_bins1.end(), vec_bins2.begin(), vec_bins2.end(), std::back_inserter(intersection));

  // Calculate the number of shared bins considering the bin frequencies
  int B12 = 0;
  for (int bin : intersection)
    B12 += std::min(bins1[bin], bins2[bin]);

  // Calculate the fraction of shared bins
  double fraction_shared = static_cast<double>(2 * B12) / (spec1.size() + spec2.size());


  return fraction_shared > min_shared_ion_fraction;
}

map<double, vector<int>> NeighborSeq::createMassPositonMap(const vector<AASequence>& candidates)
{
  map<double, vector<int>> mass_position_map;

  //Iterate through the vector of AASequence objects 
  for (size_t i = 0; i < candidates.size(); ++i)
  {
    // Calculate the mono-isotopic mass of the sequence
    double mass = candidates[i].getMonoWeight();
    // Insert the mass and the position into the map
    mass_position_map[mass].push_back(i);
  }
  return mass_position_map;
}

// Method to find neighbor peptides in a given FASTA file
vector<int> NeighborSeq::findNeighborPeptides(const AASequence& peptides,
                                                                const vector<AASequence>& neighbor_candidate,
                                                                const vector<int>& candidate_position,
                                                                const double& min_shared_ion_fraction,
                                                                const double& mz_bin_size)

{
    vector<int> result_entries;
    String rel_sequence = peptides.toString();
    MSSpectrum spec = generateSpectrum(rel_sequence);
    for (size_t i = 0; i < candidate_position.size(); i++)
    {
          String neighbor_sequence = neighbor_candidate[candidate_position[i]].toString();
          //cout << neighbor_sequence << endl;
          // Check if the sequence contains an 'X'
              if (neighbor_sequence.find('X') == String::npos)
              {     
                  MSSpectrum neighbor_spec = generateSpectrum(neighbor_sequence);

                  // Compare the spectra and add to results if they are similar
                      if (compareSpectra(spec, neighbor_spec, min_shared_ion_fraction, mz_bin_size))
                      {      
                          result_entries.push_back(candidate_position[i]);
                      }
               }
    }
return result_entries;
}

