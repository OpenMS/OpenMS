import argparse
import pyopenms as pms
from common import addDataProcessing, writeParamsIfRequested, updateDefaults


def align(in_files, out_files, trafo_out_files, reference_index,
        reference_file, params):

    in_types = set(pms.FileHandler.getType(in_) for in_ in in_files)


    if in_types <= set((pms.Type.MZML, pms.Type.MZXML, pms.Type.MZDATA)):
        align_features = False
    elif in_types == set((pms.Type.FEATUREXML,)):
        align_features = True
    else:
        parser.error("different kinds of input files")

    algo = pms.MapAlignmentAlgorithmPoseClustering()
    algo.setReference(reference_index, reference_file)

    model_params = params.copy("model:", True)
    model_type   = model_params.getValue("type").toString()

    plog = pms.ProgressLogger()
    plog.setLogType(pms.LogType.CMD)

    alignment_param = params.copy("algorithm:", True)

    algo.setParameters(alignment_param)

    transformations = []

    in_types = set(pms.FileHandler.getType(in_file) for in_file in in_files)
    if reference_file is not None:
        reference_index = in_files.index(reference_file)
    else:
        if reference_index in (None, 0):
            sizes = []
            if align_features:
                fh = pms.FeatureXMLFile()
                plog.startProgress(0, len(in_files), "Determine Reference map")
                for i, in_f in enumerate(in_files):
                    sizes.append((fh.loadSize(in_f), i))
                    plog.setProgress(i)
            else:
                fh = pms.MzMLFile()
                mse = pms.MSExperiment()
                plog.startProgress(0, len(in_files), "Determine Reference map")
                for i, in_f in enumerate(in_files):
                    fh.load(in_f, mse)
                    mse.updateRanges(1)
                    sizes.append((mse.getSize(), i))
                    plog.setProgress(i)
            plog.endProgress()
            __, reference_index = max(sizes)

    assert reference_index >=0 and reference_index < len(in_files)

    ref_file = in_files[reference_index]

    if align_features:
        f_fxml = pms.FeatureXMLFile()
        options = f_fxml.getOptions()
        fm = pms.FeatureMap()
        f_fxml_tmp = pms.FeatureXMLFile()








    in_maps = []
    if in_types <= set((pms.Type.MZML, pms.Type.MZXML, pms.Type.MZDATA)):
        fh = pms.FileHandler()
        pl.startProgress(0, len(in_files), "loading input files")
        for i, in_file in enumerate(in_files):
            pl.setProgress(i)
            pm = pms.MSExperiment()
            fh.loadExperiment(in_file, pm)
            in_maps.append(pm)
        pl.endProgress()
        algo.alignPeakMaps(in_maps, transformations)
        if model_type != "none":
            algo.fitModel(model_type, model_params, transformations)
        pms.MapAlignmentAlgorithmPoseClustering.transformPeakMaps(in_maps, transformations)
        pl.startProgress(0, len(out_files), "writing output files")
        for i, out_file in enumerate(out_files):
            pl.setProgress(i)
            in_map = addDataProcessing(in_maps[i], params, pms.ProcessingAction.ALIGNMENT)
            fh.storeExperiment(out_file, in_map)
        pl.endProgress()

    elif in_types == set((pms.Type.FEATUREXML,)):
        fh = pms.FeatureXMLFile()
        pl.startProgress(0, len(in_files), "loading input files")
        for i, in_file in enumerate(in_files):
            pl.setProgress(i)
            pm = pms.FeatureMap()
            fh.load(in_file, pm)
            in_maps.append(pm)
        pl.endProgress()
        algo.alignFeatureMaps(in_maps, transformations)
        if model_type != "none":
            algo.fitModel(model_type, model_params, transformations)
        pms.MapAlignmentAlgorithmPoseClustering.transformFeatureMaps(in_maps, transformations)
        pl.startProgress(0, len(out_files), "writing output files")
        for i, out_file in enumerate(out_files):
            pl.setProgress(i)
            in_map = addDataProcessing(in_maps[i], params, pms.ProcessingAction.ALIGNMENT)
            fh.store(out_file, in_map)
        pl.endProgress()

    else:
        raise Exception("can not handle input file format")

    if trafo_out_files:
        for name, trafo in zip(trafo_out_files, transformations):
            pms.TransformationXMLFile().store(name, trafo)


def getModelDefaults(default_model):
    params = pms.Param()
    params.setValue("type", pms.DataValue(default_model), "Type of model")
    model_types = [ "linear", "b_spline", "interpolated"]
    if default_model not in model_types:
        model_types.insert(0, default_model)
    params.setValidStrings("type", model_types)

    model_params = pms.Param()

    pms.TransformationModelLinear.getDefaultParameters(model_params)
    params.insert("linear:", model_params)
    params.setSectionDescription("linear", "Parameters for 'linear' model")

    pms.TransformationModelBSpline.getDefaultParameters(model_params)
    params.insert("b_spline:", model_params)
    params.setSectionDescription("b_spline", "Parameters for 'b_spline' model")

    pms.TransformationModelInterpolated.getDefaultParameters(model_params)
    entry = model_params.getEntry("interpolation_type")
    interpolation_types = entry.valid_strings
    if "polynomial" in interpolation_types:
        interpolation_types.remove("polynomial")
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

    run_mode = args.in_ is not None and args.out is not None\
                and (args.ini is not None or args.dict_ini is not None)
    write_mode = args.write_ini is not None or args.write_dict_ini is not None
    ok = run_mode or write_mode
    if not ok:
        parser.error("either specify -in, -(trafo_)out and -(dict)ini for running "
                     "the peakpicker\nor -write(dict)ini for creating std "
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
