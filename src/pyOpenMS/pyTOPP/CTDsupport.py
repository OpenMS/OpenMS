import CTDopts
import sys
from CTDopts.CTDopts import CTDModel, parse_cl_directives
import pyopenms as pms

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

