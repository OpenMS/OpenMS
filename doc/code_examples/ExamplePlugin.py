#!/usr/bin/env python
import argparse
import os
import sys

from pyopenms import Param
from pyopenms import ParamXMLFile

'''
The tool needs basic parameter parsing. The following parameters are required to be supported:
-write_ini <filename> : create the ini file for the tool
-ini <filename> : the ini file to load
-in <filename> : the input file
-out <filename> : the output file
-no_progress : flag to disable output to command line
'''
def main():
    parser = argparse.ArgumentParser(description="Test tool for prototyping TOPPView-plugins")
    parser.add_argument("-write_ini", help="Writes ini to specified path and exits.")

    args = parser.parse_args()

    ini_path = args.write_ini

    # create ini at path
    if ini_path is not None:
        #Create the default parameters
        param = Param()
        # this will create the param structure that is mandatory for all plugins
        param.initPluginParam("ExamplePlugin", "0.0.1")
        # the valid input fileformats have to be added like this
        param.setValidStrings("ExamplePlugin:1:in", [b"*.mzML"])

        # additional parameters can be added like this
        param.setValue("ExamplePlugin:1:Number", 0, "This is an additional numeric parameter")
        param.setValue("ExamplePlugin:1:Text", "example 1", "This is an additional text parameter")
        param.setValidStrings("ExamplePlugin:1:Text", [b"option 1", b"option 2", b"option 3"])
        param.setValue("ExamplePlugin:1:Required", "", "This is an additional required parameter", [b"required"])

        #Write them to the given filepath
        file = ParamXMLFile()
        file.store(ini_path, param)

        exit()

    
if __name__ == "__main__":
    main()


