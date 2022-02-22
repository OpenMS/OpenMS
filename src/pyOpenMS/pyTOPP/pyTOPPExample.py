import CTDopts
import sys
from CTDopts.CTDopts import CTDModel, parse_cl_directives
import pyopenms as pms
import logging
from common import addDataProcessing
from CTDsupport import *


def main():

    # register command line arguments
    model = CTDModel(
       name='NameOfThePyTOPPTool',  # required
       version='1.0',  # required
       description='This is an example tool how to write pyTOPP tools compatible with the OpenMS workflow ecosystem.',
       manual='RTF',
       docurl='http://dummy.url/docurl.html',
       category='Example',
       executableName='exampletool',
       executablePath='/path/to/exec/exampletool-1.0/exampletool'
    )

    # Register in / out etc. with CTDModel
    model.add(
        'input',
        required=True,
        type='input-file',
        is_list=False,
        file_formats=['mzML'],  # filename restrictions
        description='Input file'
    )

    model.add(
        'output',
        required=True,
        type='output-file',
        is_list=False,
        file_formats=['mzML'],  # filename restrictions
        description='Output file'
    )

    defaults = pms.PeakPickerHiRes().getDefaults()

    # expose algorithm parameters in command line options
    addParamToCTDopts(defaults, model)

    # parse command line
    # if -write_ini is provided, store model in CTD file, exit with error code 0
    # if -ini is provided, load CTD file into defaults Param object and return new model with paraneters set as defaults
    arg_dict, openms_params = parseCTDCommandLine(sys.argv, model, defaults);

    # data processing
    fh = pms.MzMLFile()
    fh.setLogType(pms.LogType.CMD)
    input_map = pms.MSExperiment()

    fh.load(arg_dict["input"], input_map)

    pp = pms.PeakPickerHiRes()
    pp.setParameters(openms_params)
    out_map = pms.MSExperiment()
    pp.pickExperiment(input_map, out_map)

    out_map = addDataProcessing(out_map, openms_params, pms.DataProcessing.ProcessingAction.PEAK_PICKING)
    fh = pms.FileHandler()
    fh.storeExperiment(arg_dict["output"], out_map)

if __name__ == "__main__":
    main()

