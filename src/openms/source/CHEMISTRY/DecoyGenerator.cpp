#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CHEMISTRY/ProteaseDigestion.h>
#include <OpenMS/CHEMISTRY/DecoyGenerator.h>

#include <random>

using namespace OpenMS;

// static
AASequence DecoyGenerator::reverseProtein(const AASequence& aas)
{
    String s = aas.toUnmodifiedString();
    std::reverse(s.begin(), s.end());
    return AASequence::fromString(s);
}

// static
AASequence DecoyGenerator::reversePeptides(const AASequence& aas, const String& protease)
{
    std::vector<AASequence> peptides;
    ProteaseDigestion ed;
    ed.setMissedCleavages(0); // important as we want to reverse between all cutting sites
    ed.setEnzyme(protease);
    ed.setSpecificity(EnzymaticDigestion::SPEC_NONE);
    ed.digest(aas, peptides);    
    String pseudo_reversed;
    for (const auto & aas : peptides)
    {
        std::string s = aas.toUnmodifiedString();
        auto last = --s.end(); // don't reverse enzymatic cutting site
        std::reverse(s.begin(), last);
        pseudo_reversed += s;
    }
    return AASequence::fromString(pseudo_reversed);
}

AASequence DecoyGenerator::shufflePeptide(
        const AASequence& protein,
        const String& protease,
        const int max_attempts,
        int seed)
{
    std::random_device r;
    std::mt19937 rng(seed);

    std::vector<AASequence> peptides;
    ProteaseDigestion ed;
    ed.setMissedCleavages(0); // important as we want to reverse between all cutting sites
    ed.setEnzyme(protease);
    ed.setSpecificity(EnzymaticDigestion::SPEC_NONE);
    ed.digest(protein, peptides);    
    String protein_shuffled;
    for (const auto & aas : peptides)
    {
        std::string peptide_string = aas.toUnmodifiedString();

        String peptide_string_shuffled = peptide_string;
        auto last = --peptide_string_shuffled.end();
        double lowest_identity(1.0);
        String lowest_identity_string(peptide_string_shuffled);
        for (int i = 0; i < max_attempts; ++i) // try to find sequence with identity lower than threshold
        {
          std::shuffle(std::begin(peptide_string_shuffled), last, rng);
          double identity = SequenceIdentity_(peptide_string_shuffled, peptide_string);
          if (identity < lowest_identity)
          {
              lowest_identity = identity;
              lowest_identity_string = peptide_string_shuffled;
          }
        }

        protein_shuffled += lowest_identity_string;
    }
    return AASequence::fromString(protein_shuffled);
}

// static
double DecoyGenerator::SequenceIdentity_(const String& decoy, const String& target)
{
    int match = 0;
    for (Size i = 0; i < target.size(); ++i)
    {
        if (target[i] == decoy[i]) { ++match; }
    }
    double identity = (double) match / target.size();
    return identity;
}    
