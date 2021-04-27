import CTDopts
import sys
import os
from CTDopts.CTDopts import CTDModel, parse_cl_directives
import pyopenms as pms
import tempfile


# code related to CTD support
def addParamToCTDopts(defaults, model):
    keys = defaults.keys()
    for key in keys:
        ctd_tags = defaults.getTags(key)
        value = defaults.getValue(key)
        desc = defaults.getDescription(key)

        ctd_required = False
        if "required" in ctd_tags:
            ctd_required = True

        ctd_type = None
        ctd_type_str = ''
        ctd_list = False

        if isinstance(value, int):
            ctd_type = int
            ctd_type_str = 'int'
        elif isinstance(value, float):
            ctd_type = float
            ctd_type_str = 'float'
        elif isinstance(value, str): 
            ctd_type = str
            ctd_type_str = 'str'
        elif isinstance(value, list):
            ctd_list = True
            # we can't determine the type based on an element so we need this helper function
            value_type = defaults.getValueType(key)
            if value_type == pms.ValueType.STRING_LIST:
                ctd_type=str
                ctd_type_str = 'list of str'
            elif value_type == pms.ValueType.DOUBLE_LIST:
                ctd_type=float
                ctd_type_str = 'list of float'
            elif value_type == pms.ValueType.INT_LIST:
                ctd_type=int
                ctd_type_str = 'list of int'

        print('Adding parameter: {0} with value: {1} and description: "{2}".'.format(key, value, desc))
        print('        required: {0} \t tags: {1} \t type: {2}.'.format(ctd_required, ctd_tags, ctd_type_str))

        model.add(
            key.decode(),
            required=ctd_required,
            type=ctd_type,
            default=value,
            is_list=ctd_list,
            description=desc)


def parseCTDCommandLine(argv, model, openms_param):
    # Configure CTDOpt to use OpenMS style on the command line.
    directives = parse_cl_directives(argv, input_ctd='ini', write_tool_ctd='write_ini', prefix='-')

    if directives["write_tool_ctd"] is not None: # triggered if -write_ini was provided on CML
        # if called with -write_ini write CTD
        model.write_ctd(directives["write_tool_ctd"])
        exit(0)
    elif directives["input_ctd"] is not None: # read ctd/ini file
        model = CTDModel(from_file=directives["input_ctd"])
#        print(model.get_defaults())
           
        param = pms.Param()
        fh = pms.ParamXMLFile()
        fh.load(directives["input_ctd"], param)
        openms_param.update(param, True)
        return model.get_defaults(), openms_param
        
    else: # only command line options provided
        temp = tempfile.NamedTemporaryFile(suffix='ini') # makes sure we get a writable file
        tmp_name = temp.name
        temp.close() # removes the file

        model.write_ctd(tmp_name)
        param = pms.Param()
        fh = pms.ParamXMLFile()
        fh.load(tmp_name, param)
        openms_param.update(param)
        os.remove(tmp_name)
        return model.parse_cl_args(argv), openms_param
