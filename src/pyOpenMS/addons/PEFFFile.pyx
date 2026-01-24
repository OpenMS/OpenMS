




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
