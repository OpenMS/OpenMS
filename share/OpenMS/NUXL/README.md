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

## Example

See the default `nuxl_presets.json` in the share folder for all presets for RNA and DNA.