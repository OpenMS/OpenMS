import argparse
import pyopenms as pms
from common import addDataProcessing, writeParamsIfRequested, updateDefaults
from  collections import Counter


def id_mapper(in_file, id_file, out_file, params):

    in_type = pms.FileHandler.getType(in_file)

    if in_types == pms.Type.CONSENSUSXML:
        link_features = 0
    elif in_types == pms.Type.FEATUREXML:
        link_features = 1
    elif in_types == pms.Type.MZQ:
        link_features = 2


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

    print args



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

        id_mapper(arg.in_, args.id_, args.out, defaults)



if __name__ == "__main__":
    main()
