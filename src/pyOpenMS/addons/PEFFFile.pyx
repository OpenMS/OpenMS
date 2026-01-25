




def get_modifications_dict(self):
    """Return modifications as a dictionary of lists suitable for DataFrame creation.

    Returns:
        dict: Dictionary with keys 'position', 'accession', 'name', 'evidence'
              where each value is a list of the corresponding field values.

    Example:
        >>> entry = PEFFEntry()
        >>> # ... load entry ...
        >>> data = entry.get_modifications_dict()
        >>> import pandas as pd
        >>> df = pd.DataFrame(data)
    """
    data = {'position': [], 'accession': [], 'name': [], 'evidence': []}
    for mod in self.modifications:
        data['position'].append(mod.position)
        data['accession'].append(mod.accession)
        data['name'].append(mod.name)
        data['evidence'].append(mod.evidence)
    return data


def get_variants_dict(self):
    """Return simple variants as a dictionary of lists suitable for DataFrame creation.

    Returns:
        dict: Dictionary with keys 'position', 'variant_aa', 'sources'
              where each value is a list of the corresponding field values.

    Example:
        >>> entry = PEFFEntry()
        >>> # ... load entry ...
        >>> data = entry.get_variants_dict()
        >>> import pandas as pd
        >>> df = pd.DataFrame(data)
    """
    data = {'position': [], 'variant_aa': [], 'sources': []}
    for var in self.simple_variants:
        data['position'].append(var.position)
        data['variant_aa'].append(chr(var.variant_aa) if var.variant_aa else '')
        data['sources'].append(var.sources)
    return data


def get_processed_regions_dict(self):
    """Return processed regions as a dictionary of lists suitable for DataFrame creation.

    Returns:
        dict: Dictionary with keys 'start_position', 'end_position', 'type', 'name', 'description'
              where each value is a list of the corresponding field values.

    Example:
        >>> entry = PEFFEntry()
        >>> # ... load entry ...
        >>> data = entry.get_processed_regions_dict()
        >>> import pandas as pd
        >>> df = pd.DataFrame(data)
    """
    data = {'start_position': [], 'end_position': [], 'type': [], 'name': [], 'description': []}
    for reg in self.processed_regions:
        data['start_position'].append(reg.start_position)
        data['end_position'].append(reg.end_position)
        data['type'].append(reg.type)
        data['name'].append(reg.name)
        data['description'].append(reg.description)
    return data


def getVariantSequences(self):
    """Get all variant sequences (each simple variant applied individually).

    Each variant sequence has one amino acid substitution applied.
    Does not combine variants.

    Returns:
        list: List of tuples (variant_description, AASequence)

    Example:
        >>> entry = PEFFEntry()
        >>> # ... load entry with variants ...
        >>> for desc, seq in entry.getVariantSequences():
        ...     print(f"{desc}: {seq.getMonoWeight():.2f} Da")
    """
    from pyopenms import AASequence
    result = []
    base_seq = self.sequence

    for var in self.simple_variants:
        if var.position == 0 or var.position > len(base_seq):
            continue

        # Apply variant: replace amino acid at position (1-based)
        idx = var.position - 1
        variant_char = chr(var.variant_aa) if var.variant_aa else ''
        if not variant_char:
            continue

        variant_seq_str = base_seq[:idx] + variant_char + base_seq[idx + 1:]

        # Create description
        original_aa = base_seq[idx] if idx < len(base_seq) else '?'
        desc = f"{original_aa}{var.position}{variant_char}"
        if var.sources:
            desc += f" ({var.sources})"

        try:
            variant_seq = AASequence.fromString(variant_seq_str)
            result.append((desc, variant_seq))
        except Exception:
            # Skip if sequence can't be parsed (invalid amino acid)
            pass

    return result


def digestWithVariants(self, digestor, min_length=6, max_length=40,
                       include_reference=True, include_variants=True,
                       include_modifications=True):
    """Digest protein and enumerate PEFF variant/modification combinations.

    Args:
        digestor: ProteaseDigestion object configured with enzyme
        min_length: Minimum peptide length (default: 6)
        max_length: Maximum peptide length (default: 40, 0 = no limit)
        include_reference: Include unmodified reference peptides (default: True)
        include_variants: Include PEFF variant peptides (default: True)
        include_modifications: Include PEFF modification combinations (default: True)

    Returns:
        list: List of tuples (description, AASequence)

    Example:
        >>> from pyopenms import ProteaseDigestion, PEFFEntry
        >>> digestor = ProteaseDigestion()
        >>> digestor.setEnzyme("Trypsin")
        >>> entry = PEFFEntry()
        >>> # ... load entry ...
        >>> for desc, seq in entry.digestWithVariants(digestor):
        ...     print(f"{desc}: {seq.toString()}")
    """
    from pyopenms import AASequence

    result = []
    base_seq = self.sequence

    # Digest the protein
    protein = AASequence.fromString(base_seq)
    peptide_ranges = []
    digestor.digest(protein, peptide_ranges, min_length, max_length if max_length > 0 else 1000000)

    for start, end in peptide_ranges:
        peptide_str = base_seq[start:end]

        # Collect variants and modifications in this peptide range
        local_variants = []
        local_mods = []

        if include_variants:
            for var in self.simple_variants:
                if start < var.position <= end:
                    local_variants.append(var)

        if include_modifications:
            for mod in self.modifications:
                if mod.position > 0 and start < mod.position <= end:
                    local_mods.append(mod)

        # Generate combinations
        def enumerate_combinations(peptide_str, variants, mods, start_pos):
            """Enumerate all variant/mod combinations."""
            combinations = []

            # Reference peptide
            if include_reference or (not variants and not mods):
                combinations.append((peptide_str, ""))

            # Single variants
            for var in variants:
                idx = var.position - start_pos - 1
                if 0 <= idx < len(peptide_str):
                    variant_char = chr(var.variant_aa) if var.variant_aa else ''
                    if variant_char:
                        var_pep = peptide_str[:idx] + variant_char + peptide_str[idx + 1:]
                        orig_aa = peptide_str[idx]
                        desc = f"{orig_aa}{var.position}{variant_char}"
                        combinations.append((var_pep, desc))

            return combinations

        # Generate peptide combinations
        combos = enumerate_combinations(peptide_str, local_variants, local_mods, start)

        for pep_str, desc in combos:
            try:
                pep_seq = AASequence.fromString(pep_str)
                result.append((desc, pep_seq))
            except Exception:
                pass

    return result


def generatePeptides(self, digestor, fixed_mods=None, variable_mods=None,
                     max_variable_mods_per_peptide=2, min_length=6, max_length=40,
                     include_reference=True, include_peff_variants=True,
                     include_peff_modifications=True):
    """Generate peptides with PEFF annotations and sample handling modifications.

    Full workflow combining:
    1. Enzymatic digestion
    2. PEFF variant enumeration
    3. PEFF modification enumeration
    4. Fixed sample modifications (e.g., Carbamidomethyl)
    5. Variable sample modifications (e.g., Oxidation)

    Args:
        digestor: ProteaseDigestion object configured with enzyme
        fixed_mods: List of fixed modification names (e.g., ["Carbamidomethyl (C)"])
        variable_mods: List of variable modification names (e.g., ["Oxidation (M)"])
        max_variable_mods_per_peptide: Max variable mods per peptide (default: 2)
        min_length: Minimum peptide length (default: 6)
        max_length: Maximum peptide length (default: 40)
        include_reference: Include unmodified reference peptides (default: True)
        include_peff_variants: Include PEFF variant peptides (default: True)
        include_peff_modifications: Include PEFF modifications (default: True)

    Returns:
        list: List of tuples (description, AASequence)

    Example:
        >>> from pyopenms import ProteaseDigestion, PEFFEntry
        >>> digestor = ProteaseDigestion()
        >>> digestor.setEnzyme("Trypsin")
        >>> entry = PEFFEntry()
        >>> # ... load entry ...
        >>> peptides = entry.generatePeptides(
        ...     digestor,
        ...     fixed_mods=["Carbamidomethyl (C)"],
        ...     variable_mods=["Oxidation (M)"],
        ...     max_variable_mods_per_peptide=2
        ... )
        >>> for desc, seq in peptides:
        ...     print(f"{desc}: {seq.toString()}")
    """
    from pyopenms import AASequence, ModifiedPeptideGenerator

    if fixed_mods is None:
        fixed_mods = []
    if variable_mods is None:
        variable_mods = []

    # Get base peptides with PEFF variants
    base_peptides = self.digestWithVariants(
        digestor, min_length, max_length,
        include_reference, include_peff_variants, include_peff_modifications
    )

    result = []

    # Prepare modification maps if needed
    fixed_mods_map = None
    var_mods_map = None

    if fixed_mods:
        fixed_mods_map = ModifiedPeptideGenerator.getModifications(fixed_mods)
    if variable_mods:
        var_mods_map = ModifiedPeptideGenerator.getModifications(variable_mods)

    for desc, peptide in base_peptides:
        # Apply fixed modifications
        if fixed_mods_map is not None:
            ModifiedPeptideGenerator.applyFixedModifications(fixed_mods_map, peptide)

        if var_mods_map is None or not variable_mods:
            result.append((desc, peptide))
        else:
            # Generate variable modification combinations
            var_peptides = []
            ModifiedPeptideGenerator.applyVariableModifications(
                var_mods_map, peptide, max_variable_mods_per_peptide, var_peptides, True
            )

            for var_pep in var_peptides:
                full_desc = desc
                if var_pep.isModified() and var_pep.toString() != peptide.toString():
                    if full_desc:
                        full_desc += ";+varmod"
                    else:
                        full_desc = "+varmod"
                result.append((full_desc, var_pep))

    return result
