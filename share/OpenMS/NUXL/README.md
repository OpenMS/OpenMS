# NuXL Custom Presets

This directory contains custom presets for the OpenNuXL tool. These presets define the nucleotide configurations and fragment adducts used in the cross-linking analysis.

## JSON Format

The JSON file should contain a dictionary where each key is a preset name and each value is a dictionary with the following keys:

- `target_nucleotides`: List of nucleotides and their empirical formulas (e.g., `"U=C9H13N2O9P"`)
- `mapping`: List of nucleotide mappings (e.g., `"U->U"`)
- `can_cross_link`: String of nucleotides that can form cross-links (e.g., `"U"` for RNA or `"T"` for DNA)
- `modifications`: List of modifications that can be applied to nucleotides (e.g., `"U:"`, `"U:-H2O"`)
- `fragment_adducts`: List of fragment adducts that can be generated from nucleotides (e.g., `"U:C9H10N2O5;U-H3PO4"`)
- `standalone` (optional): String of target nucleotides that only form precursor adducts on their own (e.g., `"s"` in `RNA-SDA`, where `s` is the cross-linker bound to the peptide without nucleotide). They are never combined with other nucleotides and are not restricted by `NuXL:sequence`. They need a `target_nucleotides` entry, `modifications` and `fragment_adducts` (their adducts are only scored with these), and must be listed in `can_cross_link`; no `mapping` entry is needed. OpenNuXL stops with an error otherwise. Hits whose adduct consists only of standalone nucleotides are peptide modifications, not nucleic acid cross-links: they are reported with `NuXL:isXL` = 0 and go through the peptide FDR. With `NuXL:only_xl`, they are not searched.
- `localization` (optional): How cross-links are localized (`NuXL:localization`). `"fragments"` (default): at the residue the shifted fragment ions point to. `"NHS"` (e.g., `RNA-SDA`): for cross-linkers attached through an NHS ester, by its chemistry among the residues the fragment ions can't tell apart: the free protein N-terminus, then lysine, then serine, threonine or tyrosine with a histidine within four residues (Mädler et al., J. Mass Spectrom. 2009, 44, 694–706), skipping residues with a modification other than a pure isotope substitution (e.g., `Acetyl (K)`, `TMT6plex (K)`; a pure isotope substitution such as `Label:13C(6) (K)` changes the isotopic composition, not the elemental formula) and, if the enzyme cleaves after lysine (e.g., trypsin), a C-terminal lysine unless the peptide ends the protein (the enzyme doesn't cleave after a modified lysine); otherwise the cross-link stays unlocalized. How the site was chosen is reported in `NuXL:localization_method`, the residues the fragment ions support in `NuXL:localization_range`.

## Example

See the default `nuxl_presets.json` in the share folder for all presets for RNA and DNA.