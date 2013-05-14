import argparse
import pyopenms as pms
import logging
from common import addDataProcessing, writeParamsIfRequested, updateDefaults


def run_peak_picker(input_map, params, out_path):

    spec0 = input_map[0]
    if pms.PeakTypeEstimator().estimateType(spec0) == \
        pms.                               SpectrumSettings.SpectrumType.PEAKS:
        logging.warn("input peak map does not look like profile data")

    if any(not s.isSorted() for s in input_map):
        raise Exception("Not all spectra are sorted according to m/z")

    pp = pms.PeakPickerHiRes()
    pp.setParameters(params)
    out_map = pms.MSExperiment()
    pp.pickExperiment(input_map, out_map)

    out_map = addDataProcessing(out_map, params, pms.ProcessingAction.PEAK_PICKING)
    fh = pms.FileHandler()
    fh.storeExperiment(out_path, out_map)


def main():

    parser = argparse.ArgumentParser(description="PeakPickerHiRes")
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

    args = parser.parse_args()

    run_mode = args.in_ is not None and args.out is not None\
                and (args.ini is not None or args.dict_ini is not None)
    write_mode = args.write_ini is not None or args.write_dict_ini is not None
    ok = run_mode or write_mode
    if not ok:
        parser.error("either specify -in, -out and -(dict)ini for running "
                     "the peakpicker\nor -write(dict)ini for creating std "
                     "ini file")

    defaults = pms.PeakPickerHiRes().getDefaults()

    write_requested = writeParamsIfRequested(args, defaults)

    if not write_requested:
        updateDefaults(args, defaults)

        fh = pms.MzXMLFile()
        fh.setLogType(pms.LogType.CMD)
        input_map = pms.MSExperiment()
        fh.load(args.in_, input_map)

        run_peak_picker(input_map, defaults, args.out)

if __name__ == "__main__":
    main()
