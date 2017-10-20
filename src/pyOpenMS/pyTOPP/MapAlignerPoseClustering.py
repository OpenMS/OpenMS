import argparse
import pyopenms as pms
from common import addDataProcessing, writeParamsIfRequested, updateDefaults


def align(in_files, out_files, out_trafos, reference_index,
          reference_file, params):

    in_types = set(pms.FileHandler.getType(in_) for in_ in in_files)

    if in_types <= set((pms.Type.MZML, pms.Type.MZXML, pms.Type.MZDATA)):
        align_features = False
    elif in_types == set((pms.Type.FEATUREXML,)):
        align_features = True
    else:
        raise Exception("different kinds of input files")

    algorithm = pms.MapAlignmentAlgorithmPoseClustering()
    alignment_params = params.copy("algorithm:", True)
    algorithm.setParameters(alignment_params)
    algorithm.setLogType(pms.LogType.CMD)

    plog = pms.ProgressLogger()
    plog.setLogType(pms.LogType.CMD)

    if reference_file:
        file_ = reference_file
    elif reference_index > 0:
        file_ = in_files[reference_index-1]
    else:
        sizes = []
        if align_features:
            fh = pms.FeatureXMLFile()
            plog.startProgress(0, len(in_files), "Determine Reference map")
            for i, in_f in enumerate(in_files):
                sizes.append((fh.loadSize(in_f), in_f))
                plog.setProgress(i)
        else:
            fh = pms.MzMLFile()
            mse = pms.MSExperiment()
            plog.startProgress(0, len(in_files), "Determine Reference map")
            for i, in_f in enumerate(in_files):
                fh.load(in_f, mse)
                mse.updateRanges()
                sizes.append((mse.getSize(), in_f))
                plog.setProgress(i)
        plog.endProgress()
        __, file_ = max(sizes)

    f_fmxl = pms.FeatureXMLFile()
    if not out_files:
        options = f_fmxl.getOptions()
        options.setLoadConvexHull(False)
        options.setLoadSubordinates(False)
        f_fmxl.setOptions(options)

    if align_features:
        map_ref = pms.FeatureMap()
        f_fxml_tmp = pms.FeatureXMLFile()
        options = f_fmxl.getOptions()
        options.setLoadConvexHull(False)
        options.setLoadSubordinates(False)
        f_fxml_tmp.setOptions(options)
        f_fxml_tmp.load(file_, map_ref)
        algorithm.setReference(map_ref)
    else:
        map_ref = pms.MSExperiment()
        pms.MzMLFile().load(file_, map_ref)
        algorithm.setReference(map_ref)

    plog.startProgress(0, len(in_files), "Align input maps")
    for i, in_file in enumerate(in_files):
        trafo = pms.TransformationDescription()
        if align_features:
            map_ = pms.FeatureMap()
            f_fxml_tmp = pms.FeatureXMLFile()
            f_fxml_tmp.setOptions(f_fmxl.getOptions())
            f_fxml_tmp.load(in_file, map_)
            if in_file == file_:
                trafo.fitModel("identity")
            else:
                algorithm.align(map_, trafo)
            if out_files:
                pms.MapAlignmentTransformer.transformRetentionTimes(map_, trafo)
                addDataProcessing(map_, params, pms.ProcessingAction.ALIGNMENT)
                f_fxml_tmp.store(out_files[i], map_)
        else:
            map_ = pms.MSExperiment()
            pms.MzMLFile().load(in_file, map_)
            if in_file == file_:
                trafo.fitModel("identity")
            else:
                algorithm.align(map_, trafo)
            if out_files:
                pms.MapAlignmentTransformer.transformRetentionTimes(map_, trafo)
                addDataProcessing(map_, params, pms.ProcessingAction.ALIGNMENT)
                pms.MzMLFile().store(out_files[i], map_)
        if out_trafos:
            pms.TransformationXMLFile().store(out_trafos[i], trafo)

        plog.setProgress(i+1)

    plog.endProgress()



def getModelDefaults(default_model):
    params = pms.Param()
    params.setValue("type", default_model, "Type of model")
    model_types = [ "linear", "interpolated"]
    if default_model not in model_types:
        model_types.insert(0, default_model)
    params.setValidStrings("type", model_types)

    model_params = pms.Param()

    pms.TransformationModelLinear.getDefaultParameters(model_params)
    params.insert("linear:", model_params)
    params.setSectionDescription("linear", "Parameters for 'linear' model")

    pms.TransformationModelInterpolated.getDefaultParameters(model_params)
    entry = model_params.getEntry("interpolation_type")
    interpolation_types = entry.valid_strings
    model_params.setValidStrings("interpolation_type", interpolation_types)

    params.insert("interpolated:", model_params)
    params.setSectionDescription("interpolated", "Parameters for 'interpolated' model")
    return params


def getDefaultParameters():
    model_param = getModelDefaults("linear")
    algo_param = pms.MapAlignmentAlgorithmPoseClustering().getParameters()
    default = pms.Param()
    default.insert("model:", model_param)
    default.insert("algorithm:", algo_param)
    return default


def main():

    parser = argparse.ArgumentParser(description="PeakPickerHiRes")
    parser.add_argument("-in",
                        action="append",
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
                        action="append",
                        type=str,
                        metavar="output_file",
                        )

    parser.add_argument("-trafo_out",
                        action="append",
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

    parser.add_argument("-reference:file",
                        action="store",
                        type=str,
                        metavar="reference_file",
                        dest="reference_file",
                        )
    parser.add_argument("-reference:index",
                        action="store",
                        type=int,
                        metavar="reference_index",
                        dest="reference_index",
                        )

    args = parser.parse_args()

    def collect(args):
        return [f.strip() for arg in args or [] for f in arg.split(",")]

    in_files = collect(args.in_)
    out_files = collect(args.out)
    trafo_out_files = collect(args.trafo_out)

    run_mode = (in_files and (out_files or trafo_out_files))\
                and (args.ini is not None or args.dict_ini is not None)

    write_mode = args.write_ini is not None or args.write_dict_ini is not None
    ok = run_mode or write_mode
    if not ok:
        parser.error("either specify -in, -(trafo_)out and -(dict)ini for running "
                     "the map aligner\nor -write(dict)ini for creating std "
                     "ini file")

    defaults = getDefaultParameters()
    write_requested = writeParamsIfRequested(args, defaults)

    if not write_requested:
        updateDefaults(args, defaults)

        if not out_files and not trafo_out_files:
            parser.error("need -out or -trafo_out files")

        if out_files and len(out_files) != len(in_files):
            parser.error("need as many -out files as -in files")
        if trafo_out_files and len(trafo_out_files) != len(in_files):
            parser.error("need as many -trafo_out files as -in files")

        if args.reference_index is not None and args.reference_file is not None:
            parser.error("can only handle either reference:index or reference:file")

        if args.reference_index is not None:
            if args.reference_index <0 or args.reference_index >= len(in_files):
                parser.error("reference:index invalid")
        if args.reference_file is not None:
            if args.reference_file not in in_files:
                parser.error("reference_file not in input files")


        align(in_files, out_files, trafo_out_files, args.reference_index or 0,
                args.reference_file or "", defaults)



if __name__ == "__main__":
    main()
