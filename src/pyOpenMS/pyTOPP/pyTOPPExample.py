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

    # Configure CTDOpt to use OpenMS style on the command line.
    directives = parse_cl_directives(sys.argv, input_ctd='ini', write_tool_ctd='write_ini', prefix='-')

    if directives["write_tool_ctd"] is not None: # triggered if -write_ini was provided on CML
        # if called with -write_ini write CTD
        model.write_ctd(directives["write_tool_ctd"])
        exit(0)
    elif directives["input_ctd"] is not None: # read ctd/ini file
        param = pms.Param()
        fh = pms.ParamXMLFile()
        fh.load(directives["input_ctd"], param)
        defaults.update(param)
     
    # data processing
    fh = pms.MzMLFile()
    fh.setLogType(pms.LogType.CMD)
    input_map = pms.MSExperiment()
    fh.load(model.input, input_map)

    pp = pms.PeakPickerHiRes()
    pp.setParameters(params)
    out_map = pms.MSExperiment()
    pp.pickExperiment(input_map, out_map)

    out_map = addDataProcessing(out_map, params, pms.ProcessingAction.PEAK_PICKING)
    fh = pms.FileHandler()
    fh.storeExperiment(model.output, out_map)

if __name__ == "__main__":
    main()
