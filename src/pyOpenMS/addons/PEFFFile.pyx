




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
