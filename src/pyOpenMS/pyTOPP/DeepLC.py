import sys
from CTDopts.CTDopts import CTDModel
from CTDsupport import addParamToCTDopts, parseCTDCommandLine
import pyopenms as pms
from io import StringIO


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

    input_csv = idxml_to_csv(arg_dict["input"])


def idxml_to_csv(idxml_file: str) -> str:
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
    return csv_str.getvalue()


if __name__ == "__main__":
    main()
