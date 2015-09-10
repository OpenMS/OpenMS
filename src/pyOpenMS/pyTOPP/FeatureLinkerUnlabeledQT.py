import argparse
import pyopenms as pms
from common import addDataProcessing, writeParamsIfRequested, updateDefaults
from  collections import Counter


def link(in_files, out_file, keep_subelements, params):

    in_types = set(pms.FileHandler.getType(in_) for in_ in in_files)

    if in_types == set((pms.Type.CONSENSUSXML,)):
        link_features = False
    elif in_types == set((pms.Type.FEATUREXML,)):
        link_features = True
    else:
        raise Exception("different kinds of input files")

    algorithm_parameters = params.copy("algorithm:", True)
    algorithm = pms.FeatureGroupingAlgorithmQT()
    algorithm.setParameters(algorithm_parameters)

    out_map = pms.ConsensusMap()
    fds = out_map.getFileDescriptions()
    if link_features:
        f = pms.FeatureXMLFile()
        maps = []
        for i, in_file in enumerate(in_files):
            map_ = pms.FeatureMap()
            f.load(in_file, map_)

            # set filedescriptions
            fd = fds.get(i, pms.FileDescription())
            fd.filename = in_file
            fd.size = map_.size()
            fd.unique_id = map_.getUniqueId()
            fds[i] = fd
            maps.append(map_)
        out_map.setFileDescriptions(fds)
        algorithm.group(maps, out_map)
    else:
        f = pms.ConsensusXMLFile()
        maps = []
        for i, in_file in enumerate(in_files):
            map_ = pms.ConsensusMap()
            f.load(in_file, map_)
            maps.append(map_)
        algorithm.group(maps, out_map)

        if not keep_subelements:
            for i in range(len(in_files)):
                # set filedescriptions
                fd = fds.get(i, pms.FileDescription())
                fd.filename = in_files[i]
                fd.size = maps[i].size()
                fd.unique_id = maps[i].getUniqueId()
                fds[i] = fd
            out_map.setFileDescriptions(fds)
        else:
            algorithm.transferSubelements(maps, out_map)

    out_map.setUniqueIds()
    addDataProcessing(out_map, params, pms.ProcessingAction.FEATURE_GROUPING)

    pms.ConsensusXMLFile().store(out_file, out_map)

    sizes = []
    for feat in out_map:
        sizes.append(feat.size())

    c = Counter(sizes)
    print "Number of consensus features:"
    for size, count in c.most_common():
        print "   of size %2d : %6d" % (size, count)
    print "        total : %6d" % out_map.size()


def main():

    parser = argparse.ArgumentParser(description="FeatureLinkerUnlabeledQT")
    parser.add_argument("-in",
                        action="append",
                        type=str,
                        dest="in_",
                        metavar="input_files",
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

    parser.add_argument("-keep_subelements",
                        action="store_true",
                        )


    args = parser.parse_args()


    def collect(args):
        return [f.strip() for arg in args or [] for f in arg.split(",")]

    in_files = collect(args.in_)

    run_mode = (in_files and args.out) \
                and (args.ini is not None or args.dict_ini is not None)

    write_mode = args.write_ini is not None or args.write_dict_ini is not None
    ok = run_mode or write_mode
    if not ok:
        parser.error("either specify -in, -out and -(dict)ini for running "
                     "the feature linker\nor -write(dict)ini for creating std "
                     "ini file")

    defaults = pms.FeatureGroupingAlgorithmQT().getParameters()
    write_requested = writeParamsIfRequested(args, defaults)

    if not write_requested:
        updateDefaults(args, defaults)


        link(in_files, args.out, args.keep_subelements, defaults)



if __name__ == "__main__":
    main()
