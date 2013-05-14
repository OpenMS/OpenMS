from special_autowrap_conversionproviders import *
from autowrap.ConversionProvider import special_converters

def register_converters():

    special_converters.append(OpenMSStringConverter())
    special_converters.append(StdVectorStringConverter())
    special_converters.append(StdSetStringConverter())
    special_converters.append(OpenMSIntListConverter())
    special_converters.append(OpenMSStringListConverter())
    special_converters.append(OpenMSDoubleListConverter())
    special_converters.append(OpenMSMapConverter())
    special_converters.append(CVTermMapConverter())
    special_converters.append(OpenMSDataValue())
    special_converters.append(OpenMSDPosition2())
    special_converters.append(OpenMSDPosition2Vector())
