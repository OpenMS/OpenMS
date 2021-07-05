import sys
from CTDopts.CTDopts import CTDModel
from CTDsupport import addParamToCTDopts, parseCTDCommandLine
import pyopenms as pms
from deeplc import DeepLC
import pandas as pd


def main():
    # register command line arguments
    model = CTDModel(
        name="DeepLC",
        version="0.1.29",
        description="DeepLC: Retention time prediction for (modified) peptides using Deep Learning.",
        docurl="",
        category="",
        executableName="",
        executablePath=""
    )
    # Register in / out etc. with CTDModel
    model.add(
        "input",
        required=True,
        type="input-file",
        is_list=False,
        file_formats=["idxml"],
        description="Input file"
    )
    model.add(
        "output",
        required=True,
        type="output-file",
        is_list=False,
        file_formats=["idxml"],
        description="Output file"
    )

    # TODO: add flag to create calibrants to model

    defaults = {}
    addParamToCTDopts(defaults, model)

    # parse command line
    # if -write_ini is provided, store model in CTD file, exit with error code 0
    # if -ini is provided, load CTD file into defaults Param object and return new model with paraneters set as defaults
    arg_dict, openms_params = parseCTDCommandLine(sys.argv, model, defaults)

    infile = arg_dict["input"]
    outfile = arg_dict["output"]

    # load search engine results
    pep_ids = []
    prot_ids = []
    pms.IdXMLFile().load(infile, prot_ids, pep_ids)

    # TODO: if calibrants flag is set use calibrants

    # annotate predictions
    prot_ids, pep_ids = annotate_predictions(prot_ids, pep_ids, outfile, predict(peptide_ids_to_dataframe(pep_ids)))

    #calibrants = create_calibration_data(pep_ids)
    #prot_ids, pep_ids = annotate_predictions(prot_ids, pep_ids, outfile, predict(peptide_ids_to_dataframe(pep_ids), calibrants))


def peptide_ids_to_dataframe(pep_ids: list) -> pd.DataFrame:
    """Parse a given list of peptide identification to a pandas DataFrame compatible with DeepLC.

    See https://github.com/compomics/DeepLC#input-files for more information.

    Args:
        pep_ids: List containing PeptideIdentification.

    Returns:
        pandas DataFrame with three columns: seq, modifications and tr.
        The returned pandas DataFrame can be used directly by DeepLC.
    """
    columns = ["seq", "modifications", "tr"]
    sequences = []
    modifications = []
    tr = []
    for pep_id in pep_ids:
        rt = pep_id.getRT()
        for hit in pep_id.getHits():
            seq = hit.getSequence()
            sequences.append(seq.toUnmodifiedString())
            tr.append(rt)
            hit_mods = []
            for pos in range(0, seq.size()):
                residue = seq.getResidue(pos)
                if residue.isModified():
                    hit_mods.append("|".join([str(pos + 1), residue.getModificationName()]))
            modifications.append("|".join(hit_mods))
    data = {
        "seq": sequences,
        "modifications": modifications,
        "tr": tr
    }
    return pd.DataFrame(data, columns=columns)


def create_calibration_data(pep_ids: list) -> pd.DataFrame:
    """Create pandas DataFrame to be used by DeepLC for calibration.

    Args:
        pep_ids: List containing PeptideIdentification
    """
    # check if target_decoy information present
    has_target_decoy = False
    for pep_id in pep_ids:
        for hit in pep_id.getHits():
            if hit.metaValueExists("target_decoy"):
                has_target_decoy = True
    if has_target_decoy:
        # annotate q-value
        fdr = pms.FalseDiscoveryRate()
        fdr.apply(pep_ids)

        # filter by 1% PSM FDR (q < 0.01)
        idfilter = pms.IDFilter()
        idfilter.filterHitsByScore(pep_ids, 0.01)
        idfilter.removeDecoyHits(pep_ids)

    return peptide_ids_to_dataframe(pep_ids)


def annotate_predictions(prot_ids: list, pep_ids: list, predictions: list) -> list:
    """Annotates predictions and returns proteins and peptides.

    Adds a custom meta value containing the prediction made by DeepLC.

    Args:
        prot_ids:       the protein identifications
        pep_ids:        the peptide identifications (PSMs)
        predictions:    A list predictions to be inserted into the resulting idXML file.
    """

    preds_iter = iter(predictions)
    # Insert prediction for each hit of a peptide id as meta value
    for pep_id in pep_ids:
        new_hits = []
        for hit in pep_id.getHits():
            try:
                hit.setMetaValue("prediction", str(next(preds_iter)))
            except StopIteration:
                raise SystemExit("Error: Number of predictions and peptide hits does not match.")
            new_hits.append(hit)
        pep_id.setHits(new_hits)
    return [prot_ids, pep_is]


def predict(peptides: pd.DataFrame, calibration: pd.DataFrame = None) -> list:
    """
    Make predications based on a given peptide DataFrame and an optional
    calibration DataFrame using DeepLC.
    """
    dlc = DeepLC()
    if calibration is None:
        return dlc.make_preds(seq_df=peptides, calibrate=False)
    else:
        dlc.calibrate_preds(seq_df=calibration)
        return dlc.make_preds(seq_df=peptides)


if __name__ == "__main__":
    main()
