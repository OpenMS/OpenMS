import argparse
import pyopenms as pms
from common import addDataProcessing, writeParamsIfRequested, updateDefaults


def id_mapper(in_file, id_file, out_file, params, use_centroid_rt,
        use_centroid_mz, use_subelements ):

    in_type = pms.FileHandler.getType(in_file)

    protein_ids = []
    peptide_ids = []

    pms.IdXMLFile().load(id_file, protein_ids, peptide_ids)

    mapper = pms.IDMapper()
    mapper.setParameters(params)

    if in_type == pms.Type.CONSENSUSXML:
        file_ = pms.ConsensusXMLFile()
        map_ = pms.ConsensusMap()
        file_.load(in_file, map_)
        mapper.annotate(map_, peptide_ids, protein_ids, use_subelements)
        addDataProcessing(map_, params, pms.DataProcessing.ProcessingAction.IDENTIFICATION_MAPPING)
        file_.store(out_file, map_)

    elif in_type == pms.Type.FEATUREXML:
        file_ = pms.FeatureXMLFile()
        map_ = pms.FeatureMap()
        file_.load(in_file, map_)
        mapper.annotate(map_, peptide_ids, protein_ids, use_centroid_rt,
                use_centroid_mz)
        addDataProcessing(map_, params, pms.DataProcessing.ProcessingAction.IDENTIFICATION_MAPPING)
        file_.store(out_file, map_)

    else:
        raise Exception("invalid input file format")


def main():

    parser = argparse.ArgumentParser(description="IDMapper")

    parser.add_argument("-id",
                        action="store",
                        type=str,
                        dest="id_",
                        metavar="id_file",
                        )

    parser.add_argument("-in",
                        action="store",
                        type=str,
                        dest="in_",
                        metavar="input_file",
                        )

    parser.add_argument("-out",
                        action="store",
                        type=str,
                        metavar="output_file",
                        )

    parser.add_argument("-ini",
                        action="store",
                        type=str,
                        metavar="ini_file",
                        )

    parser.add_argument("-dict_ini",
                        action="store",
                        type=str,
                        metavar="python_dict_ini_file",
                        )

    parser.add_argument("-write_ini",
                        action="store",
                        type=str,
                        metavar="ini_file",
                        )

    parser.add_argument("-write_dict_ini",
                        action="store",
                        type=str,
                        metavar="python_dict_ini_file",
                        )

    parser.add_argument("-mz_tolerance",
                        action="store",
                        type=float)

    parser.add_argument("-rt_tolerance",
                        action="store",
                        type=float)

    parser.add_argument("-mz_measure",
                        action="store",
                        type=str)

    parser.add_argument("-mz_reference",
                        action="store",
                        type=str)

    parser.add_argument("-ignore_charge",
                        action="store_true")

    parser.add_argument("-feature:use_centroid_rt",
                        action="store_true")

    parser.add_argument("-feature:use_centroid_mz",
                        action="store_true")

    parser.add_argument("-consensusfeature:use_subelements",
                        action="store_true")

    args = parser.parse_args()



    run_mode = (args.in_ and args.id_ and args.out) \
                and (args.ini is not None or args.dict_ini is not None)

    write_mode = args.write_ini is not None or args.write_dict_ini is not None
    ok = run_mode or write_mode
    if not ok:
        parser.error("either specify -id, -in, -out and -(dict)ini for running "
                     "the id mapper\nor -write(dict)ini for creating std "
                     "ini file")

    defaults = pms.IDMapper().getParameters()
    write_requested = writeParamsIfRequested(args, defaults)

    if not write_requested:
        updateDefaults(args, defaults)

        dd = defaults.asDict()
        dd["ignore_charge"] = "true" if args.ignore_charge else "false"
        for att in ["mz_tolerance", "rt_tolerance", "mz_measure",
                "mz_reference"]:
            value = getattr(args, att)
            if value is not None:
                dd[att] = value

        defaults.update(dd)
        feature_use_centroid_rt = getattr(args, "feature:use_centroid_rt")
        feature_use_centroid_mz = getattr(args, "feature:use_centroid_mz")
        consenususfeature_use_subelements = getattr(args,
                "consensusfeature:use_subelements")

        id_mapper(args.in_, args.id_, args.out, defaults,
                feature_use_centroid_rt,
                feature_use_centroid_mz,
                consenususfeature_use_subelements
                )



if __name__ == "__main__":
    main()
