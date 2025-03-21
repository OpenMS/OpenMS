// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow, Philipp Wang $
// $Authors: Chris Bielow, Philipp Wang $
// --------------------------------------------------------------------------
#include <OpenMS/ANALYSIS/ID/NeighborSeq.h>
#include <OpenMS/CHEMISTRY/ResidueDB.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/KERNEL/MSSpectrum.h>

#include <OpenMS/MATH/MathFunctions.h>

using namespace OpenMS;
using namespace std;


NeighborSeq::NeighborSeq(std::vector<AASequence>&& digested_relevant_peptides)
  : digested_relevant_peptides_(std::move(digested_relevant_peptides)),
    neighbor_stats_(digested_relevant_peptides_.size(), 0)
{
  Param params;
  params.setValue("add_b_ions", "true");
  params.setValue("add_y_ions", "true");
  params.setValue("add_first_prefix_ion", "true"); // do not skip b1 ion
  spec_gen_.setParameters(params);

  x_residue_ = ResidueDB::getInstance()->getResidue('X');

  // Index peptide masses for fast lookup
  mass_position_map_ = createMassLookup_();
}

// Function to generate the theoretical spectrum for a given peptide sequence
MSSpectrum NeighborSeq::generateSpectrum(const AASequence& peptide_sequence)
{
  MSSpectrum spectrum;
  spec_gen_.getSpectrum(spectrum, peptide_sequence, 1, 1);
  return spectrum;
}

int NeighborSeq::computeSharedIonCount(const MSSpectrum& spec1, const MSSpectrum& spec2, const double& mz_bin_size)
{
  // compute shared b/y ions in two sorted ranges
  auto setIntersectionCount = [mz_bin_size](auto first1, auto last1, auto first2, auto last2) -> Size
  {
    Size count {0};
    while (first1 != last1 && first2 != last2)
    {
      auto val1 = int(first1->getMZ() / mz_bin_size);
      auto val2 = int(first2->getMZ() / mz_bin_size);
      if (val1 < val2) ++first1;
      else
      {
        if (val1 == val2)
        {
          ++first1;
          ++count;
        }
        ++first2;
      }
    }
    return count;
  };

  auto shared_ions = setIntersectionCount(spec1.begin(), spec1.end(), spec2.begin(), spec2.end());

  return shared_ions;
}

// Function to compare two spectra and determine if they are similar
bool NeighborSeq::isNeighborSpectrum(const MSSpectrum& spec1, const MSSpectrum& spec2, const double min_shared_ion_fraction, const double mz_bin_size)
{
  // Calculate the number of shared bins considering the bin frequencies
  int B12 = computeSharedIonCount(spec1,  spec2, mz_bin_size);

  // Calculate the fraction of shared bins
  double fraction_shared = (2.0 * B12) / (spec1.size() + spec2.size());

  return fraction_shared > min_shared_ion_fraction;
}

//Finds candidate positions based on a given mono-isotopic weight and mass tolerance.
auto NeighborSeq::findCandidatePositions_(const double mono_weight, double mass_tolerance, const bool mass_tolerance_pc_ppm)
{
  // Calculate the lower and upper bounds for the mass tolerance range
  assert(mass_tolerance >= 0);
  if (mass_tolerance_pc_ppm)
  {
    mass_tolerance = Math::ppmToMass(mono_weight, mass_tolerance);
  }

  // Find the lower bound iterator in the map
  auto lower = mass_position_map_.lower_bound(mono_weight - mass_tolerance);

  // Find the upper bound iterator in the map
  auto upper = mass_position_map_.upper_bound(mono_weight + mass_tolerance);

  return make_pair(lower, upper);
}

// Method to find neighbor peptides in a given FASTA file
bool NeighborSeq::isNeighborPeptide(const AASequence& peptide,
                                    const double mass_tolerance_pc,
                                    const bool mass_tolerance_pc_ppm,
                                    const double min_shared_ion_fraction,
                                    const double mz_bin_size)

{
  auto [from, to] = findCandidatePositions_(peptide.getMonoWeight(), mass_tolerance_pc, mass_tolerance_pc_ppm);
  if (from == to) return false;

  bool found = false;
  MSSpectrum spec = generateSpectrum(peptide);
  for (auto it_rel_pep = from; it_rel_pep != to; ++it_rel_pep)
  {
    for (int pep_index : it_rel_pep->second)
    {
      MSSpectrum neighbor_spec = generateSpectrum(digested_relevant_peptides_[pep_index]);
      if (isNeighborSpectrum(spec, neighbor_spec, min_shared_ion_fraction, mz_bin_size))
      { 
        //std::cout << digested_relevant_peptides_[pep_index] << " has neighbor " << peptide << '\n';
        neighbor_stats_[pep_index]++;
        found = true;
      }
    }
  }
  return found;
}

map<double, vector<int>> NeighborSeq::createMassLookup_()
{
  // Map to store the mass and corresponding positions
  map<double, vector<int>> mass_position_map;

  int skipped{0};
  // Iterate through the vector of AASequence objects
  for (size_t i = 0; i < digested_relevant_peptides_.size(); ++i)
  {
    if (digested_relevant_peptides_[i].has(*x_residue_))
    {
      neighbor_stats_[i] = -1; // mark as not findable
      skipped++;
      continue;
    }
    // Calculate the mono-isotopic mass of the sequence
    double mass = digested_relevant_peptides_[i].getMonoWeight();

    // Insert the mass and the position into the map
    mass_position_map[mass].push_back(i);
  }
  OPENMS_LOG_WARN << "Skipped " << skipped << "/" << digested_relevant_peptides_.size()
                  << " peptides with unknown('X') amino acids." << endl;
  return mass_position_map;
}

NeighborSeq::NeighborStats NeighborSeq::getNeighborStats() const
{
  NeighborStats stats;
  for (int count : neighbor_stats_)
  {
    if (count == -1)
      stats.unfindable_peptides++;
    else if (count == 0)
      stats.findable_no_neighbors++;
    else if (count == 1)
      stats.findable_one_neighbor++;
    else
      stats.findable_multiple_neighbors++;
  }
  return stats;
}




