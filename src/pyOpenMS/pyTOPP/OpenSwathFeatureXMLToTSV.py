
import pyopenms

def convert_to_row(first, targ, run_id, keys, filename):
    peptide_ref = first.getMetaValue("PeptideRef")
    pep = targ.getPeptideByRef(peptide_ref)
    full_peptide_name = "NA"
    if (pep.metaValueExists("full_peptide_name")):
      full_peptide_name = pep.getMetaValue("full_peptide_name")

    decoy = "0"
    peptidetransitions = [t for t in targ.getTransitions() if t.getPeptideRef() == peptide_ref]
    if len(peptidetransitions) > 0:
        if peptidetransitions[0].getDecoyTransitionType() == pyopenms.DecoyTransitionType().DECOY:
            decoy = "1"
        elif peptidetransitions[0].getDecoyTransitionType() == pyopenms.DecoyTransitionType().TARGET:
            decoy = "0"

    protein_name = "NA"
    if len(pep.protein_refs) > 0:
      protein_name = pep.protein_refs[0]

    row = [
        first.getMetaValue("PeptideRef"),
        run_id,
        filename,
        first.getRT(),
        first.getMetaValue("PrecursorMZ"),
        first.getUniqueId(),
        pep.sequence,
        full_peptide_name,
        pep.getChargeState(),
        first.getMetaValue("PrecursorMZ"),
        first.getIntensity(),
        protein_name,
        decoy
    ]

    for k in keys:
        row.append(first.getMetaValue(k))

    return row

def get_header(features):
    keys = []
    features[0].getKeys(keys)
    header = [
        "transition_group_id",
        "run_id",
        "filename",
        "RT",
        "id",
        "Sequence" ,
        "FullPeptideName",
        "Charge",
        "m/z",
        "Intensity",
        "ProteinName",
        "decoy"]
    header.extend(keys)
    return header
