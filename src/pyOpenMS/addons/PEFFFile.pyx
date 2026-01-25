




def get_modifications_dict(self):
    """Return modifications as dict for DataFrame creation."""
    data = {'position': [], 'accession': [], 'name': [], 'evidence': []}
    for mod in self.modifications:
        data['position'].append(mod.position)
        data['accession'].append(mod.accession)
        data['name'].append(mod.name)
        data['evidence'].append(mod.evidence)
    return data


def get_variants_dict(self):
    """Return simple variants as dict for DataFrame creation."""
    data = {'position': [], 'variant_aa': [], 'sources': []}
    for var in self.simple_variants:
        data['position'].append(var.position)
        data['variant_aa'].append(chr(var.variant_aa) if var.variant_aa else '')
        data['sources'].append(var.sources)
    return data


def get_processed_regions_dict(self):
    """Return processed regions as dict for DataFrame creation."""
    data = {'start_position': [], 'end_position': [], 'type': [], 'name': [], 'description': []}
    for reg in self.processed_regions:
        data['start_position'].append(reg.start_position)
        data['end_position'].append(reg.end_position)
        data['type'].append(reg.type)
        data['name'].append(reg.name)
        data['description'].append(reg.description)
    return data


def getVariantSequences(self):
    """Get variant sequences (each simple variant applied individually).

    Returns list of (description, AASequence) tuples.
    """
    from pyopenms import AASequence
    result = []
    base_seq = self.sequence

    for var in self.simple_variants:
        if var.position == 0 or var.position > len(base_seq):
            continue
        idx = var.position - 1
        variant_char = chr(var.variant_aa) if var.variant_aa else ''
        if not variant_char:
            continue

        variant_seq_str = base_seq[:idx] + variant_char + base_seq[idx + 1:]
        orig_aa = base_seq[idx] if idx < len(base_seq) else '?'
        desc = f"{orig_aa}{var.position}{variant_char}"
        if var.sources:
            desc += f" ({var.sources})"

        try:
            result.append((desc, AASequence.fromString(variant_seq_str)))
        except Exception:
            pass

    return result


def digestWithVariants(self, digestor, min_length=6, max_length=40,
                       include_reference=True, include_variants=True):
    """Digest protein and apply PEFF variants.

    Returns list of (description, AASequence) tuples.
    """
    from pyopenms import AASequence

    result = []
    base_seq = self.sequence
    protein = AASequence.fromString(base_seq)
    peptide_ranges = []
    digestor.digest(protein, peptide_ranges, min_length, max_length if max_length > 0 else 1000000)

    for start, end in peptide_ranges:
        peptide_str = base_seq[start:end]

        # Add reference peptide
        if include_reference:
            try:
                result.append(("", AASequence.fromString(peptide_str)))
            except Exception:
                pass

        # Add variant peptides
        if include_variants:
            for var in self.simple_variants:
                if start < var.position <= end:
                    idx = var.position - start - 1
                    variant_char = chr(var.variant_aa) if var.variant_aa else ''
                    if variant_char and 0 <= idx < len(peptide_str):
                        var_pep = peptide_str[:idx] + variant_char + peptide_str[idx + 1:]
                        desc = f"{peptide_str[idx]}{var.position}{variant_char}"
                        try:
                            result.append((desc, AASequence.fromString(var_pep)))
                        except Exception:
                            pass

    return result


def generatePeptides(self, digestor, fixed_mods=None, variable_mods=None,
                     max_variable_mods_per_peptide=2, min_length=6, max_length=40,
                     include_reference=True, include_peff_variants=True):
    """Generate peptides with PEFF variants and sample modifications.

    Args:
        digestor: ProteaseDigestion object
        fixed_mods: Fixed mods like ["Carbamidomethyl (C)"]
        variable_mods: Variable mods like ["Oxidation (M)"]

    Returns list of (description, AASequence) tuples.
    """
    from pyopenms import ModifiedPeptideGenerator

    base_peptides = self.digestWithVariants(
        digestor, min_length, max_length, include_reference, include_peff_variants
    )

    if not fixed_mods and not variable_mods:
        return base_peptides

    result = []
    fixed_map = ModifiedPeptideGenerator.getModifications(fixed_mods) if fixed_mods else None
    var_map = ModifiedPeptideGenerator.getModifications(variable_mods) if variable_mods else None

    for desc, peptide in base_peptides:
        if fixed_map:
            ModifiedPeptideGenerator.applyFixedModifications(fixed_map, peptide)

        if not var_map:
            result.append((desc, peptide))
        else:
            var_peptides = []
            ModifiedPeptideGenerator.applyVariableModifications(
                var_map, peptide, max_variable_mods_per_peptide, var_peptides, True
            )
            for var_pep in var_peptides:
                full_desc = desc
                if var_pep.toString() != peptide.toString():
                    full_desc = f"{desc};+varmod" if desc else "+varmod"
                result.append((full_desc, var_pep))

    return result
