// Compare OpenDIA::readLibrary with TransitionPQPFile (field by field, same order). Validation only.
#include "OpenDIALibrary.h"
#include <OpenMS/ANALYSIS/OPENSWATH/TransitionPQPFile.h>
#include <iostream>
using namespace OpenMS;
int main(int argc, char** argv)
{
  for (int a = 1; a < argc; ++a)
  {
    OpenSwath::LightTargetedExperiment ref;
    TransitionPQPFile().convertPQPToTargetedExperiment(argv[a], ref);
    const auto lib = OpenDIA::readLibrary(argv[a]);
    std::size_t bad = 0;
    if (ref.compounds.size() != lib.compounds.size() || ref.transitions.size() != lib.transitions.size()) bad++;
    for (std::size_t i = 0; i < std::min(ref.compounds.size(), lib.compounds.size()); ++i)
    {
      const auto& r = ref.compounds[i]; const auto& l = lib.compounds[i];
      const bool rt_same = (r.rt == l.rt) || (r.rt != r.rt && l.rt != l.rt);
      if (r.id != l.id || !rt_same || r.drift_time != l.drift_time || r.charge != l.charge || r.sequence != l.sequence ||
          r.protein_refs != l.protein_refs) { if (bad++ < 3) std::cerr << "compound " << i << " " << r.id << "\n"; }
    }
    for (std::size_t i = 0; i < std::min(ref.transitions.size(), lib.transitions.size()); ++i)
    {
      const auto& r = ref.transitions[i]; const auto& l = lib.transitions[i];
      if (r.transition_name != l.transition_name || r.peptide_ref != l.peptide_ref || r.library_intensity != l.library_intensity ||
          r.product_mz != l.product_mz || r.precursor_mz != l.precursor_mz || r.precursor_im != l.precursor_im ||
          r.fragment_charge != l.fragment_charge || r.getDecoy() != l.getDecoy() || r.getAnnotation() != l.getAnnotation())
      { if (bad++ < 3) std::cerr << "transition " << i << " " << r.transition_name << " " << r.getAnnotation() << "/" << l.getAnnotation() << "\n"; }
    }
    const auto tr = TransitionPQPFile().getPQPCurrentIDToTraMLIDMap(argv[a], "PRECURSOR");
    if (tr != OpenDIA::readPrecursorTraMLIDs(argv[a])) bad++;
    std::cout << argv[a] << ": compounds " << lib.compounds.size() << ", transitions " << lib.transitions.size()
              << ", mismatches " << bad << "\n";
  }
}
