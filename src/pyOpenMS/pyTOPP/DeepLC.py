import sys
from CTDopts.CTDopts import CTDModel
from CTDsupport import addParamToCTDopts, parseCTDCommandLine
import pyopenms as pms
from io import StringIO
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

    defaults = {}
    addParamToCTDopts(defaults, model)

    # parse command line
    # if -write_ini is provided, store model in CTD file, exit with error code 0
    # if -ini is provided, load CTD file into defaults Param object and return new model with paraneters set as defaults
    arg_dict, openms_params = parseCTDCommandLine(sys.argv, model, defaults)

    infile = arg_dict["input"]
    outfile = arg_dict["output"]
    peptide_df = idxml_to_dataframe(infile)
    write_predictions_to_idxml(infile, outfile, predict(peptide_df))


def idxml_to_csv(idxml_file: str) -> StringIO:
    """Parse a given idxml file to a csv compatible with DeepLC.

    See https://github.com/compomics/DeepLC#input-files for more information.

    Args:
        idxml_file: The path of the idXML file to be parsed to a DeepLC csv.

    Returns:
        StringIO containing the csv string.
    """
    pep_ids = []
    pms.IdXMLFile().load(idxml_file, [], pep_ids)

    csv_str: StringIO = StringIO()
    csv_str.write("seq,modifications,tr\n")
    for pep_id in pep_ids:
        pep_rt = pep_id.getRT()
        for hit in pep_id.getHits():
            seq = hit.getSequence()
            csv_str.write("".join([seq.toUnmodifiedString(), ","]))
            # Find all modifications. Store them as pos|mod where pos is the
            # position of the residue and mod the modification at that position.
            modifications = []
            for pos in range(0, seq.size()):
                residue = seq.getResidue(pos)
                if residue.isModified():
                    modifications.append("|".join([str(pos + 1), residue.getModificationName()]))
            # separate each modifications with '|': seq,pos_1|mod_1|pos_2|mod_2|...,tr
            csv_str.write("|".join(modifications))
            csv_str.write("".join([",", str(pep_rt), "\n"]))

    print(csv_str.getvalue())
    assert len(pep_ids) == len(csv_str.getvalue().splitlines()) - 1
    return csv_str


def idxml_to_dataframe(idxml_file: str) -> pd.DataFrame:
    """Parse a given idxml file to a pandas DataFrame compatible with DeepLC.

    See https://github.com/compomics/DeepLC#input-files for more information.

    Args:
        idxml_file: The path of the idXML file to be parsed to a pandas DataFrame.

    Returns:
        pandas DataFrame with three columns: seq, modifications and tr containing
        the data read from idxml_file.
    """
    pep_ids = []
    pms.IdXMLFile().load(idxml_file, [], pep_ids)

    columns = ["seq", "modifications", "tr"]
    sequences = []
    modifications = []
    tr = []
    for pep_id in pep_ids:
        tr.append(pep_id.getRT())
        for hit in pep_id.getHits():
            seq = hit.getSequence()
            sequences.append(seq.toUnmodifiedString())
            m = []
            for pos in range(0, seq.size()):
                residue = seq.getResidue(pos)
                if residue.isModified():
                    m.append("|".join([str(pos + 1), residue.getModificationName()]))
            modifications.append("|".join(m))
    data = {
        "seq": sequences,
        "modifications": modifications,
        "tr": tr
    }
    return pd.DataFrame(data, columns=columns)


def write_predictions_to_idxml(infile: str, outfile: str, predictions: list):
    """Write a new idXML file to `outfile` with the given predictions.

    Load an existing idXML file (`infile`), copy and modify the peptide identification
    data by adding a custom meta value containing the prediction made by DeepLC
    and store the data as a new idXML file at `outfile`.

    Args:
        infile:     Path of the idXML file to be loaded.
        outfile:    Path where the resulting idXML will be stored.
        predictions:      A list predictions to be inserted into the resulting idXML file.
    """
    pep_ids = []
    prot_ids = []
    pms.IdXMLFile().load(infile, prot_ids, pep_ids)
    # Insert prediction for each peptide id as meta value
    assert len(pep_ids) == len(predictions), "The number of peptides and predictions does not match."
    for idx, pep_id in enumerate(pep_ids, 0):
        pep_id.setMetaValue("prediction", str(predictions[idx]))

    pms.IdXMLFile().store(outfile, prot_ids, pep_ids)


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
        dlc.make_preds(seq_df=peptides)


if __name__ == "__main__":
    main()
