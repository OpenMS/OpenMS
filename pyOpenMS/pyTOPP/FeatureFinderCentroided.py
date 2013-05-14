import argparse
import pyopenms as pms
from common import addDataProcessing, writeParamsIfRequested, updateDefaults


def run_featurefinder_centroided(input_path, params, seeds, out_path):

    fh = pms.MzMLFile()
    options = pms.PeakFileOptions()
    options.setMSLevels([1,1])
    fh.setOptions(options)
    input_map = pms.MSExperiment()
    fh.load(input_path, input_map)
    input_map.updateRanges()

    ff = pms.FeatureFinder()
    ff.setLogType(pms.LogType.CMD)

    features = pms.FeatureMap()
    name = pms.FeatureFinderAlgorithmPicked.getProductName()
    ff.run(name, input_map, features, params, seeds)

    features.setUniqueIds()
    addDataProcessing(features, params, pms.ProcessingAction.QUANTITATION)

    fh = pms.FeatureXMLFile()
    fh.store(out_path, features)


def main():

    parser = argparse.ArgumentParser(description="FeatureFinderCentroided")
    parser.add_argument("-in",
                        action="store",
                        type=str,
                        dest="in_",
                        metavar="input_file",
                        )

    parser.add_argument("-seeds",
                        action="store",
                        type=str,
                        metavar="seeds_file",
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

    args = parser.parse_args()

    run_mode = args.in_ is not None and args.out is not None\
                and (args.ini is not None or args.dict_ini is not None)
    write_mode = args.write_ini is not None or args.write_dict_ini is not None
    ok = run_mode or write_mode
    if not ok:
        parser.error("either specify -in, -out and -(dict)ini for running "
                     "the peakpicker\nor -write(dict)ini for creating std "
                     "ini file")

    name = pms.FeatureFinderAlgorithmPicked.getProductName()
    defaults = pms.FeatureFinder().getParameters(name)

    write_requested = writeParamsIfRequested(args, defaults)

    if not write_requested:
        updateDefaults(args, defaults)

        seeds = pms.FeatureMap()
        if args.seeds:
            fh = pms.FeatureXMLFile()
            fh.load(args.seeds, seeds)

        run_featurefinder_centroided(args.in_, defaults, seeds, args.out)

if __name__ == "__main__":
    main()
