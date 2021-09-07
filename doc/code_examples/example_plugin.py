#!/usr/bin/env python
import argparse
import os
import sys
from lxml import etree

NONMSP_ATTR = "xsi:noNamespaceSchemaLocation"
NONMSP_VAL = "https://raw.githubusercontent.com/OpenMS/OpenMS/develop/share/OpenMS/SCHEMAS/Param_1_7_0.xsd" 


def get_ini_xml():
    NS = 'http://www.w3.org/2001/XMLSchema-instance'
    
    location_attribute = '{%s}noNamespaceSchemaLocation' % NS

    root = etree.Element("PARAMETERS", 
        nsmap={
            "xsi" : NS,
        },
        attrib={
            location_attribute : NONMSP_VAL
        })

    tool = etree.SubElement(root, "NODE")
    tool.set("name", "ExamplePlugin")
    tool.set("description", "ExamplePlugin")

    item = etree.SubElement(tool, "ITEM")
    item.set("name", "version")
    item.set("value", "0.1")
    item.set("type", "string")
    item.set("description", "Version of the tool that generated this parameters file.")
    item.set("required", "false")
    item.set("advanced", "true")

    params = etree.SubElement(tool, "NODE")
    params.set("name", "1")
    params.set("description", "Instance '1' section for ExampleTool")

    input = etree.SubElement(params, "ITEM")
    input.set("name", "in")
    input.set("value", "")
    input.set("type", "input-file")
    input.set("description", "Input File")
    input.set("required", "false")
    input.set("advanced", "true")
    input.set("supported_formats", "*.mzid,*.idXML")

    output = etree.SubElement(params, "ITEM")
    output.set("name", "out")
    output.set("value", "")
    output.set("type", "output-file")
    output.set("description", "Output File")
    output.set("required", "false")
    output.set("advanced", "true")
    output.set("supported_formats", "*.featureXML,*.consensusXML,*.mzq")


    root.append(tool)

    return etree.ElementTree(root)


def write_ini(path):
    if len(path) < 0:
        raise ValueError("path cannot be empty")

    abs_path = os.path.abspath(path)
    if os.path.isdir(abs_path):
        raise ValueError("path cannot be a directory")
    
    with open(abs_path, "w+"):
        ini = get_ini_xml()
        ini.write(abs_path, pretty_print=True, xml_declaration=True, encoding='UTF-8')


'''
The tool needs basic parameter parsing. At the very least you should include the parsing of the -write_ini flag, that 
returns a XML ini.
'''
parser = argparse.ArgumentParser(description="Test tool for prototyping TOPPView-plugins")
parser.add_argument("-write_ini", help="Writes ini to specified path and exits.")

args = parser.parse_args()

ini_path = args.write_ini


# create ini at path
if ini_path is not None:
    write_ini(ini_path)
    sys.exit(0)
    