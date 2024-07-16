// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Sven Nahnsen $
// $Authors: Sven Nahnsen, Andreas Bertsch, Chris Bielow, Philipp Wang $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/CHEMISTRY/ProteaseDB.h>
#include <OpenMS/CHEMISTRY/ProteaseDigestion.h>
#include <OpenMS/ANALYSIS/OPENSWATH/MRMDecoy.h>
#include <OpenMS/CHEMISTRY/DigestionEnzyme.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/MATH/MathFunctions.h>
#include <OpenMS/DATASTRUCTURES/FASTAContainer.h>
#include <regex>
#include <OpenMS/ANALYSIS/ID/NeighborSeq.h>
#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGenerator.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <cmath>
#include <set>


using namespace OpenMS;
using namespace std;


// Function to calculate the mass of a peptide
double NeighborSeq::calculateMass(const AASequence& peptide)
{
return peptide.getMonoWeight();
}


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
bool NeighborSeq::compareSpectra(const MSSpectrum& spec1, const MSSpectrum& spec2, const double& ion_tolerance, const string& resolution)
{
// Determine the bin size based on resolution
double bin_size;
if (resolution == "high") { bin_size = 1.0005079; }
else { bin_size = 0.05; }

// Maps to hold the binned spectra
map<int, int> bins1, bins2;

// Lambda function to bin a spectrum
auto bin_spectrum = [&](const MSSpectrum& spec, map<int, int>& bins) {
    for (const auto& peak : spec)
    {
    int bin = static_cast<int>(peak.getMZ() / bin_size);
    bins[bin]++;
    }
};

// Bin both spectra
bin_spectrum(spec1, bins1);
bin_spectrum(spec2, bins2);

// Calculate the size of the spectra and the number of shared bins
int B1 = spec1.size();
int B2 = spec2.size();
int B12 = 0;

for (const auto& bin : bins1)
{
    if (bins2.find(bin.first) != bins2.end()) { B12 += min(bin.second, bins2[bin.first]); }
}

// Calculate the fraction of shared bins
double fraction_shared = static_cast<double>(2 * B12) / (B1 + B2);
// cout << fraction_shared << endl;
return fraction_shared > ion_tolerance;
}

// Method to find neighbor peptides in a given FASTA file
vector<FASTAFile::FASTAEntry> NeighborSeq::findNeighborPeptides(const AASequence& peptides,
                                                                const vector<FASTAFile::FASTAEntry>& neighbor_file,
                                                                const double& mass_tolerance,
                                                                const double& ion_tolerance,
                                                                const String& resolution)

{
vector<FASTAFile::FASTAEntry> result_entries;

// Calculate the mass and generate the spectrum for the input peptide
size_t size = peptides.size();
double mass = calculateMass(peptides);
String peptide_sequence = peptides.toString();
MSSpectrum spec = generateSpectrum(peptide_sequence);
// Iterate through the neighbor file entries
for (const auto& entry : neighbor_file)
{
    for (size_t offset = 0; offset + size <= entry.sequence.size(); ++offset)
    {
    String neighbor_sequence = entry.sequence.substr(offset, size);
    // Check if the sequence contains an 'X'
    if (neighbor_sequence.find('X') == String::npos)
    {
        double neighbor_mass = calculateMass(AASequence::fromString(neighbor_sequence));
        double mass_diff = std::abs(mass - neighbor_mass) / ((mass + neighbor_mass) / 2);
        // Check if the mass difference is within the tolerance
        if (mass_diff <= mass_tolerance)
        {
        // cout << mass_diff << endl;
        MSSpectrum neighbor_spec = generateSpectrum(neighbor_sequence);
        // Compare the spectra and add to results if they are similar
        if (compareSpectra(spec, neighbor_spec, ion_tolerance, resolution))
        {
            // Construct FASTAFile::FASTAEntry and add to results
            FASTAFile::FASTAEntry neighbor_entry(entry.identifier, entry.description, neighbor_sequence);
            result_entries.push_back(neighbor_entry);
        }
        }
    }
    }
}
return result_entries;
}

