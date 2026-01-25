




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
        except:
            # Skip if sequence can't be parsed (invalid amino acid)
            pass

    return result
