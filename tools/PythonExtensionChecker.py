#!/usr/bin/env python
# -*- coding: utf-8  -*-
"""
=========================================================================
        msproteomicstools -- Mass Spectrometry Proteomics Tools
=========================================================================

Copyright (c) 2013, ETH Zurich
For a full list of authors, refer to the file AUTHORS.

This software is released under a three-clause BSD license:
 * Redistributions of source code must retain the above copyright
   notice, this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright
   notice, this list of conditions and the following disclaimer in the
   documentation and/or other materials provided with the distribution.
 * Neither the name of any author or any participating institution
   may be used to endorse or promote products derived from this software
   without specific prior written permission.
--------------------------------------------------------------------------
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
--------------------------------------------------------------------------
$Maintainer: Hannes Roest$
$Authors: Hannes Roest$
--------------------------------------------------------------------------
"""

# https://pypi.python.org/pypi/breathe
# $ sudo pip install breathe

import glob, os, sys
import re, time
import argparse
from PythonCheckerLib import parse_pxd_file
from PythonCheckerLib import create_pxd_file_map
import breathe.parser
from breathe.parser.doxygen.compound import parse as doxygen_parse
from Cython.Compiler.Nodes import CEnumDefNode, CppClassNode, CTypeDefNode, CVarDefNode, CImportStatNode, CDefExternNode
from autowrap.PXDParser import CppClassDecl, CTypeDefDecl, MethodOrAttributeDecl, EnumDecl
import yaml

# Matching function
def handle_member_definition(mdef, pxd_class, cnt):
    """ Matches a doxygen member definition (mdef) to a Cython pxd file.

    """
    tres = TestResult()
    protection = mdef.get_prot() # DoxProtectionKind: public, protected, private, package
    kind = mdef.get_kind() # DoxMemberKind: define property event variable typedef enum function signal prototype friend dcop slot 
    if not protection in "public protected private package".split(" "):
        raise Exception("Error; something is wrong")
    if not kind in "define property event variable typedef enum function signal prototype friend dcop slot".split(" "):
        raise Exception("Error; something is wrong")
    if kind == "enum" and protection == "public":
        cnt.public_enums_total += 1
        cython_file = parse_pxd_file(pxd_class.pxd_path)
        found = False
        for klass in cython_file:
            if hasattr(klass[0], "name") and klass[0].name == mdef.get_name():
                found = True
                break

        if not found:
            tres.setPassed(False)
            tres.setMessage("TODO: Found enum in C++ but not in pxd: %s %s %s" % (mdef.kind, mdef.prot, mdef.name))
            cnt.public_enums_missing += 1
            comp_name = mdef.parent_doxy_file.compound.get_compoundname()
            file_location = mdef.parent_doxy_file.getCompoundFileLocation()
            internal_file_name = "OpenMS" + file_location.split("/include/OpenMS")[1]
            namespace = comp_name
            true_cppname = '"%s::%s"' % (comp_name, mdef.get_name())
            enumr  = "\n"
            enumr += 'cdef extern from "<%s>" namespace "%s":\n' % (internal_file_name, namespace)
            enumr += "    \n"
            enumr += '    cdef enum %s %s:\n' % (mdef.get_name(), true_cppname)
            for val in mdef.get_enumvalue():
                enumr += "        %s\n" % val.get_name()
            tres.setMessage(tres.getMessage() + enumr)
        elif len(klass[0].items) != len(mdef.get_enumvalue()):
            tres.setPassed(False)
            tres.setMessage("TODO: Found enum in C++ with %s members but in Cython there are %s members: " % (
                len(mdef.get_enumvalue()), len(klass[0].items) ) )
        else:
            tres.setPassed(True)

    elif kind == "variable" and protection == "public":
        attrnames = [a.name for a in pxd_class.attributes]
        cnt.public_variables += 1
        if not mdef.name in attrnames:
            tres.setPassed(False)
            tres.setMessage("TODO: Found attribute in C++ but not in pxd: %s %s %s" % (mdef.kind, mdef.prot, mdef.name) )
            cnt.public_variables_missing += 1
        else:
            tres.setPassed(True)

    elif kind == "function" and protection == "public":
        # Wrap of public member functions ... 
        cnt.public_methods += 1
        c_return_type = mdef.resolve_return_type()
        if mdef.name in pxd_class.methods:
            # Found match between C++ method and Python method
            py_methods = pxd_class.methods[mdef.name]
            if not isinstance(py_methods, list):
                py_methods = [py_methods]

            py_return_type = [str(d.result_type) for d in py_methods]
            if mdef.definition == mdef.name:
                # Constructor, no return type -> all is good
                assert len(c_return_type) == 0
                tres.setPassed(True)
            elif "void" in py_return_type and not "void" in c_return_type:
                tres.setPassed(False)
                tres.setMessage( "TODO: Mismatch between C++ return type (%s) and Python return type (%s) in %s %s %s:" % ( 
                     str(c_return_type), str(py_return_type), mdef.kind,  mdef.prot, mdef.name) )
            else:
                tres.setPassed(True)

        else:
            cnt.public_methods_missing += 1
            if mdef.name.find("~") != -1:
                # destructor
                cnt.public_methods_missing_nowrapping += 1
                tres.setPassed(True)
                tres.setMessage("Cannot wrap destructor")
            elif (mdef.name.find("operator") != -1 or
                  mdef.name.find("begin") != -1 or
                  mdef.name.find("end") != -1):
                cnt.public_methods_missing_nowrapping += 1
                tres.setPassed(True)
                tres.setMessage("Cannot wrap method with iterator/operator %s" % mdef.name)
            else:
                tres.setPassed(False)
                tres.setMessage(" -- TODO missing function in PXD: %s " % mdef.format_definition_for_cython())

    else:
        # It is neither public function/enum/variable
        tres.setPassed(True)

    # Return the testresult
    return tres
#
## Class for counting occurances
# 
class Counter(object):
    def __init__(self):
        self.total = 0
        self.skipped = 0
        self.skipped_could_not_parse = 0
        self.skipped_ignored = 0
        self.skipped_protected = 0
        self.skipped_no_location = 0
        self.skipped_no_sections = 0
        self.skipped_no_pxd_file = 0
        self.skipped_no_pxd_match = 0
        self.parsed = 0
        # 
        self.public_enums_total = 0
        self.public_enums_missing = 0
        # 
        self.public_methods = 0
        self.public_methods_missing = 0
        self.public_methods_missing_nowrapping = 0
        self.public_variables = 0
        self.public_variables_missing = 0

    def computed_skipped(self):
        self.skipped = self.skipped_could_not_parse +\
            self.skipped_ignored + \
            self.skipped_protected + \
            self.skipped_no_location + \
            self.skipped_no_sections + \
            self.skipped_no_pxd_file + \
            self.skipped_no_pxd_match  

    def print_skipping_reason(self):
        self.computed_skipped()
        print "Skipped files: %s" % self.skipped
        print "- Could not parse xml: %s" % self.skipped_could_not_parse
        print "- Could not parse location in xml: %s" % self.skipped_no_location
        print "- Ignored per ignore-file: %s" % self.skipped_ignored
        print "- Protected Compound: %s" % self.skipped_protected
        print "- Could not find sections in xml: %s" % self.skipped_no_sections
        print "- Could not find associated pxd file : %s" % self.skipped_no_pxd_file
        print "- Could not find matching class in pxd file : %s" % self.skipped_no_pxd_match

    def print_stats(self):
        self.computed_skipped()
        print "Total files: %s" % self.total
        print "Skipped files: %s" % self.skipped
        print "Parsed files: %s" % self.parsed
        print "Parsed public methods %s (of which were missing %s and %s were operator/destructors) " % (self.public_methods, self.public_methods_missing, self.public_methods_missing_nowrapping)
        print "  - wrapped %s " % (self.public_methods - self.public_methods_missing)
        print "  - unwrapped operators/destructors %s " % (self.public_methods_missing_nowrapping)
        print "  - unwrapped methods %s " % (self.public_methods_missing - self.public_methods_missing_nowrapping)
        print "Parsed public enums %s (of which were missing %s) " % (self.public_enums_total, self.public_enums_missing)
        print "Parsed public attributes %s (of which were missing %s) " % (self.public_variables, self.public_variables_missing)
        print "Note that this script counts each method name only once and only maps from \n"+ \
              "C++ to Python (not the other way around), thus the numbers are slightly inaccurate." 

#
## Class for an OpenMS .h file
# 
class OpenMSSourceFile(object):
    def __init__(self, fname):
        self.fname = fname

    def getMaintainer(self):
        try:
            return self._getMaintainer()
        except IOError:
            return None

    def _getMaintainer(self):
        """
        // $Maintainer: xxx $
        // $Authors: xxx $
        """
        maintainer_reg = re.compile(".*\$\s*Maintainer:([^\$]*)\$")
        with open(self.fname) as f:
            data = f.read()
            maintainer = maintainer_reg.search(data)
            if maintainer is not None:
                return maintainer.group(1).strip()
            else:
                return None


#
## Class holding Doxygen XML file and next one function declaration
# 
class DoxygenXMLFile(object):
    """The doxygen XML file
    """

    def __init__(self, fname):
        self.fname = fname
        self.parsed_file = None
        self.compound = None
        self.parsing_error = False
        self.parsing_error_message = None

    def parse_doxygen(self):
        try:
            self.parsed_file =  doxygen_parse(self.fname)
            self.compound = self.parsed_file.get_compounddef()
            return self.parsed_file
        except Exception as e:
            print "Error parsing doxygen xml file", e.message
            self.parsing_error_message = e.message
            self.parsing_error = True
            return None

    def getCompoundFileLocation(self):
        location = self.parsed_file.get_compounddef().get_location() 
        if location is None:
            return None
        return os.path.realpath(location.get_file())

    def isEmpty(self, discard_defines=False):
        compound = self.compound
        if not discard_defines:
            return len(compound.get_sectiondef()) == 0

        # check whether there is more than defines and typdefs
        empty = True
        for mdef in self.iterMemberDef():
            if not mdef.get_kind() in ["define", "typedef", "slot", "signal"]:
                # DoxMemberKind: define property event variable typedef enum function signal prototype friend dcop slot 
                empty = False
        if empty and not len(compound.get_sectiondef()) == 0:
            # contains only typedefs etc
            pass
        return empty

    def get_pxd_from_class(self, dfile, file_location, xml_output_path):
        compound = dfile.compound
        comp_name = compound.get_compoundname()
        #
        # Step 1: print includes
        # 
        includes = ""
        if len(compound.get_includes()) == 1:
            try:
                reffile = os.path.join(xml_output_path, compound.get_includes()[0].get_refid() + ".xml")
                # read includes from ref file
                dreffile = DoxygenXMLFile(reffile).parse_doxygen()
                include_compound = dreffile.get_compounddef()
            except Exception as e:
                print "Error: Could not read includes from file for compound %s with error %s" % (comp_name, e.message)
                include_compound = compound
        else:
            include_compound = compound

        for inc in include_compound.get_includes():
            val = inc.getValueOf_()
            if val.startswith("OpenMS"):
                header_file = val.split("/")[-1]
                if header_file.endswith(".h"):
                    header_file = header_file[:-2]
                # We do not import certain headers since we did not wrap them in Python
                if header_file in ["Exception", "Macros", "config", "StandardTypes"]:
                    continue
                includes += "from %s cimport *\n" % header_file
        if len(file_location.split("/include/OpenMS")) != 2:
            return None
        #
        # Step 2: class definition
        # 
        internal_file_name = "OpenMS" + file_location.split("/include/OpenMS")[1]
        # for n in compound.get_inheritancegraph().get_node():
        parent_classes = [n.getValueOf_() for n in compound.basecompoundref]
        namespace = "::".join(comp_name.split("::")[:-1])
        preferred_classname = "_".join(comp_name.split("::")[1:])
        preferred_classname = comp_name.split("::")[-1]
        cldef  = "\n"
        cldef += 'cdef extern from "<%s>" namespace "%s":\n' % (internal_file_name, namespace)
        cldef += "    \n"

        inherit_txt = ""
        true_cppname = '"%s"' % comp_name
        if len(parent_classes) > 0:
            inherit_txt = "(%s)" % ",".join(parent_classes)
            true_cppname = '' # Cython does not accept this if mixed with inheritance

        if compound.templateparamlist is None:
            cldef += '    cdef cppclass %s%s %s:\n' % (preferred_classname, inherit_txt, true_cppname)
        else:
            targs = [p.get_declname() for p in compound.templateparamlist.get_param()]
            cldef += '    cdef cppclass %s[%s]%s:\n' % (preferred_classname, ",".join(targs), inherit_txt)
        if len(parent_classes) > 0:
            cldef += '        # wrap-inherits:\n'
        for p in parent_classes:
            cldef += '        #  %s\n' % (p)

        # check if it is abstract, then do not attempt to wrap it
        if self.isAbstract():
            cldef += '        # wrap-ignore\n'
            cldef += '        # ABSTRACT class\n'

        #
        # Step 3: methods
        # 
        methods = ""
        default_ctor = False
        copy_ctor = False
        enum = ""
        imports_needed = {}
        for mdef in dfile.iterMemberDef():
            if mdef.kind == "enum" and mdef.prot == "public":
                # add enums
                enum += '\n'
                enum += 'cdef extern from "<%s>" namespace "%s":\n' % (internal_file_name, comp_name)
                enum += '    cdef enum %s "%s":\n' % (mdef.name, comp_name + "::" + mdef.name)
                enum += '        #wrap-attach:\n'
                enum += '        #    %s\n' % preferred_classname
                for val in mdef.enumvalue:
                    enum += '        %s\n' % val.get_name()

            if mdef.kind == "variable" and mdef.prot == "public":
                # print "var", mdef.name
                methods += "        %s\n" % mdef.format_definition_for_cython(False)
            elif mdef.kind == "function" and mdef.prot == "public":
                if mdef.definition == mdef.name:
                    # print "we have a constructor", mdef.name, mdef.get_argsstring()
                    if mdef.get_argsstring().strip() == "()":
                        # print "have default"
                        default_ctor = True
                        continue
                    elif mdef.get_argsstring().strip().find(mdef.name) != -1 and \
                         mdef.get_argsstring().strip().find(",") == -1:
                        # print "have copy"
                        copy_ctor = True
                        continue
                if mdef.name.find("~") != -1:
                    # print "we have a deconstructor", mdef.name
                    continue
                # res += "do member function/attribute : ", mdef.kind,  mdef.prot, mdef.name
                declaration = mdef.format_definition_for_cython(False)
                DoxygenCppFunction.compute_imports(declaration, imports_needed)
                if declaration.find("operator=(") != -1:
                    # assignment operator, cannot be overriden in Python
                    continue
                methods += "        %s nogil except +\n" % declaration

        # Add the imports we need
        res  = DoxygenCppFunction.generate_imports(imports_needed)
        res += includes
        res += cldef
        if default_ctor:
            res += "        %s() nogil except +\n" % comp_name.split("::")[-1]
        if not copy_ctor:
            res += "        %s(%s) nogil except + #wrap-ignore\n" % (comp_name.split("::")[-1], comp_name.split("::")[-1])
        else:
            res += "        %s(%s) nogil except +\n" % (comp_name.split("::")[-1], comp_name.split("::")[-1])
        res += methods
        res += enum
        res += "\n"
        return res

    def iterMemberDef(self):
        """Iterate over all members of this class.
        
        We do not care about the sections defined in the documentation here."""
        for sdef in self.compound.get_sectiondef():
            for mdef_ in sdef.get_memberdef():
                mdef = DoxygenCppFunction.generate_from_obj(mdef_)
                mdef.parent_doxy_file = self
                yield mdef

    def isAbstract(self):
        for mdef in self.iterMemberDef():
                if mdef.get_argsstring().endswith("=0"):
                    return True
        return False

class DoxygenCppFunction(object):
    """ A Cpp function definition from a doxygen file"""
    def __init__(self):
        # No other code here, below is the real init method!
        self.initialize_dgencpp()

    @staticmethod
    def generate_from_obj(mdef):
        """Attaches the functionality of this object to the given input object"""
        for k,v in DoxygenCppFunction.__dict__.iteritems():
            if callable(v) and not k == "__init__": 
                import types
                mdef.__dict__[k] = types.MethodType(v, mdef)
        mdef.initialize_dgencpp()
        return mdef

    @staticmethod
    def generate_imports(imports):
        res = ""
        res += "from Types cimport *\n"
        for k in sorted(imports.keys()):
            if k == "bool":
                res += "from libcpp cimport bool\n"
            else:
                res += "from libcpp.%s cimport %s as libcpp_%s\n" % (k,k,k)
        return res

    @staticmethod
    def compute_imports(declaration, imports):
        if declaration.find("libcpp_vector") != -1:
            imports["vector"] = 0
        if declaration.find("libcpp_pair") != -1:
            imports["pair"] = 0
        if declaration.find("libcpp_map") != -1:
            imports["map"] = 0
        if declaration.find("libcpp_set") != -1:
            imports["set"] = 0
        if declaration.find("libcpp_string") != -1:
            imports["string"] = 0
        if declaration.find("bool") != -1:
            imports["bool"] = 0

    def initialize_dgencpp(self):
        pass

    def resolve_return_type(self):
        res = []
        return self._resolve_type(self.get_type().content_)

    def _resolve_type(self, mtype):
        res = []
        for c in mtype:
            val = c.getValue()
            if hasattr(val, "getValueOf_"):
                res.append(val.getValueOf_())
            else:
                res.append(val)
        return res

    def format_definition_for_cython(self, replace_nogil=True):
        """Parse a doxygen function definition and write it in Cython"""
        c_return_type = self.resolve_return_type()

        # remove default arguments, Cython doesnt like them
        arguments = re.sub("\=[^,\)]*", "", self.get_argsstring())
        function_name = self.name

        arguments = "("
        nested = False
        for i,p in enumerate(self.get_param()): 
            ptype = self._resolve_type(p.get_type().content_)
            dname = p.declname
            # ignore default arguments etc ... Cython cannot use them 
            # p.defval.content_
            # replace python keywords in argument name: except, type, lamdba, map ...
            import keyword
            if keyword.iskeyword(dname):
                dname = dname + "_"
            if dname in dir(__builtins__):
                dname = dname + "_"
            # dname = dname.replace("except", "except_").replace("type", "type_").replace("lambda", "lambda_").replace("map", "map_")
            tojoin = "".join(ptype) + " " + dname.strip()
            if tojoin.count("std::") > 2:
                nested = True
            if i == 0:
                arguments += tojoin
            else:
                arguments += ", " + "".join(ptype) + " " + dname.strip()
                
        arguments += ")"
        # arguments = (std::vector<  Int  > column_indices, std::vector<  DoubleReal  > column_values, const  String  & name)

        if len(self.get_argsstring()) == 0:
            arguments = ""

        # remove returned references
        return_type = "".join(c_return_type)
        return_type = return_type.replace("&", "")
        cpp_def = return_type + " " + function_name + arguments
        cpp_def = cpp_def.replace("///", "#")    
        cpp_def = cpp_def.replace("//", "#")    
        if replace_nogil:
            cpp_def = cpp_def.replace(";", "nogil except +")    
            cpp_def = cpp_def.replace("const;", "nogil except +")    
        else:
            cpp_def = cpp_def.replace("const;", "")
            cpp_def = cpp_def.replace(";", "")
        # TODO handle static ...  
        cpp_def = cpp_def.replace("static", "")
        cpp_def = cpp_def.replace("MSSpectrum<>", "MSSpectrum[Peak1D]")
        cpp_def = cpp_def.replace("MSChromatogram<>", "MSChromatogram[ChromatogramPeak]")
        cpp_def = cpp_def.replace("std::vector", "libcpp_vector")
        cpp_def = cpp_def.replace("std::map", "libcpp_map")
        cpp_def = cpp_def.replace("std::pair", "libcpp_pair")
        cpp_def = cpp_def.replace("std::set", "libcpp_set")
        cpp_def = cpp_def.replace("std::string", "libcpp_string")
        cpp_def = cpp_def.replace("<", "[")
        cpp_def = cpp_def.replace(">", "]")
        cpp_def = cpp_def.replace("operator[]", "operator__[]")
        cpp_def = cpp_def.replace("operator]", "operator>")
        cpp_def = cpp_def.replace("operator[", "operator<")
        cpp_def = cpp_def.replace("operator__[]", "operator[]")
        cpp_def = cpp_def.replace("const ", "")    
        cpp_def = cpp_def.replace("[ DoubleReal ]", "[ double ]")    
        cpp_def = cpp_def.replace("[ Size ]", "[ size_t ]")    
        cpp_def = cpp_def.replace("[ Int ]", "[ int ]")    
        cpp_def = cpp_def.replace("RichPeakSpectrum", "MSSpectrum[RichPeak1D]")
        cpp_def = cpp_def.replace("RichPeakMap", "MSExperiment[RichPeak1D, ChromatogramPeak]")
        cpp_def = cpp_def.replace("FeatureMap[]", "FeatureMap[Feature]")
        cpp_def = cpp_def.replace("MSSpectrum[]", "MSSpectrum[Peak1D]")
        cpp_def = cpp_def.replace("MSExperiment[]", "MSExperiment[Peak1D, ChromatogramPeak]")
        # cpp_def = cpp_def.replace("Chromatogram", "MSChromatogram[ChromatogramPeak]")
        #
        cpp_def = cpp_def.replace("PeakSpectrum", "MSSpectrum[Peak1D]")
        cpp_def = cpp_def.replace("PeakMap", "MSExperiment[Peak1D, ChromatogramPeak]")

        # Alert the user to potential problems and comment out potential
        # dangerous things (raw pointers, iterators)
        if cpp_def.find("*") != -1 or \
           cpp_def.find("::iterator") != -1:
            cpp_def = "# POINTER # " + cpp_def
        if cpp_def.find("::") != -1:
            cpp_def = "# NAMESPACE # " + cpp_def
        if self.templateparamlist is not None:
            cpp_def = "# TEMPLATE # " + cpp_def
        if nested:
            cpp_def = "# NESTED STL # " + cpp_def

        return cpp_def.strip()


#
## Class for the ignore file
# 
class IgnoreFile(object):

    def __init__(self):
        self.data = {
            "IgnoreNames" : [],
            "IgnoreMethods" : {},
        }

    def load(self, fname):
        self.fname = fname
        self.data = yaml.load(open(self.fname) )["PyOpenMSChecker"]

    def isNameIgnored(self, name):
        return name in self.data["IgnoreNames"]

    def getIgnoredMethods(self, name):
        return self.data["IgnoreMethods"].get(name, [])
        

#
## Class for the .pxd file
# 
class PXDFileParseError(Exception):
    pass

class PXDFile(object):

    def __init__(self):
        pass

    @staticmethod
    def parse_multiple_files(pxdfiles, comp_name):

        def cimport(b, _, __):
            print "cimport", b.module_name, "as", b.as_name

        handlers = { CEnumDefNode : EnumDecl.parseTree,
                     CppClassNode : CppClassDecl.parseTree,
                     CTypeDefNode : CTypeDefDecl.parseTree,
                     CVarDefNode  : MethodOrAttributeDecl.parseTree,
                     CImportStatNode  : cimport,
                     }

        found = False
        # Go through all files and all classes in those files, trying to find
        # the class whose C/C++ name matches the current compound name
        for pxdfile in pxdfiles:
            cython_file = parse_pxd_file(pxdfile)
            for klass in cython_file:
                if hasattr(klass[0], "cname"):
                    if klass[0].cname == comp_name:
                        found = True
                if found: break
            if found: break

        if not found: 
            error_str = "Could not find a match for class %s in file %s" % (comp_name, pxdfile)
            raise PXDFileParseError(error_str)

        # Check if we really have a class, then initialize it
        if isinstance(klass[0], CppClassNode):
            cl = CppClassDecl.parseTree(klass[0], klass[1], klass[2])
        else: 
            print "Something is wrong, not a class"
            raise PXDFileParseError("wrong")

        cl.pxdfile = pxdfile
        for klass in cython_file:
            handler = handlers.get(type(klass[0]))
            res = handler(klass[0], klass[1], klass[2])
            if res.annotations.has_key("wrap-attach"):
                if res.annotations["wrap-attach"] == cl.name:
                    ## attach this to the above class
                    cl.methods[res.name] = res

        return cl

class TestResult:
    """ A Result from a single test which either passed or failed.

    In addition, it contains information about the reason why the test failed,
    who is the maintainer and a unique testname.
    """

    def __init__(self, passed = None, message=None, log_level=None, name=None):
        self.message = message
        self.passed = passed
        self.name = name
        self.log_level = log_level
        self.maintainer = None
        if log_level is None:
            self.log_level = 0

    def setMessage(self, message_):
        self.message = message_

    def getMessage(self):
        return self.message

    def setPassed(self, passed):
        self.passed = passed

    def isPassed(self):
        return self.passed

    def setMaintainer(self, maintainer):
        self.maintainer = maintainer

    def getMaintainer(self):
        if self.maintainer is None:
            return "Nobody"
        return self.maintainer

class TestResultHandler:
    """ A Container for all test results.

    """

    def __init__(self):
        self._list = []

    def append(self, l):
        self._list.append(l)

    def __iter__(self):
        for l in self._list:
            yield l

    def to_cdash_xml(self, template_path, output_path):

        xml_output = []
        # load template head (everything up to "<Testing>")
        # -> this assumes a specific format of the xml
        with open(template_path) as f:
            for line in f:
                if line.strip() == "<Testing>":
                    break
                xml_output.append(line)

        """
          # load template head
          $template = file($ctestReportingPath."/Test.xml");
          $newTestFile = array();
          foreach ($template as $line)
          {
            array_push($newTestFile, $line);
            if (trim($line) == "<Testing>")
            {
              break;
            }
          }
          """

        # Start writing the xml
        xml_output.append("  <Testing>\n")
        xml_output.append("  <StartDateTime>%s</StartDateTime>\n" % (time.strftime('%b %d %H:%M') ) )
        xml_output.append("  <StartTestTime>%s</StartTestTime>\n" % (time.time()) )
        xml_output.append("  <TestList>\n")
        for classtestresults in self:
            for tres in classtestresults:
                xml_output.append("    <Test>%s</Test>\n" % tres.name) 
        xml_output.append("  </TestList>\n")
        
        for classtestresults in self:
            for tres in classtestresults:
                status = ""
                if tres.isPassed():
                    status = "passed"
                else:
                    status = "failed"

                xml_output.append(" " * 2 + '<Test Status="%s">\n' % status)
                xml_output.append(" " * 4 + '<Name>%s</Name>\n' % tres.name)
                xml_output.append(" " * 4 + '<Path> ./tools/ </Path>\n' )
                xml_output.append(" " * 4 + '<FullName>%s</FullName>\n' % tres.name)
                xml_output.append(" " * 4 + '<FullCommandLine>python PythonExtensionChecker.py %s</FullCommandLine>\n' % tres.name)
                xml_output.append(" " * 4 + '<Results>')
                xml_output.append("""
      <NamedMeasurement type="numeric/double" name="Execution Time"><Value>0.001</Value></NamedMeasurement>
      <NamedMeasurement type="text/string" name="Completion Status"><Value>Completed</Value></NamedMeasurement>
      <NamedMeasurement type="text/string" name="Maintainer"><Value>%s</Value></NamedMeasurement>
      <NamedMeasurement type="text/string" name="Command Line"><Value>python PythonExtensionChecker.py</Value></NamedMeasurement>\n""" % (
          tres.getMaintainer()
      ) )
                xml_output.append(" " * 6 + '<Measurement>\n')
                xml_output.append(" " * 8 + '<Value>\n')
                if not tres.getMessage() is None:
                    xml_output.append(" " * 10 +  tres.getMessage() + "\n")
                xml_output.append(" " * 8 + '</Value>\n')
                xml_output.append(" " * 6 + '</Measurement>\n')

                xml_output.append(" " * 4 + '</Results>\n')
                xml_output.append(" " * 2 + '</Test>\n')


        xml_output.append("<EndDateTime>%s</EndDateTime>\n" % (time.strftime('%b %d %H:%M') ) ) 
        xml_output.append("<EndTestTime>%s</EndTestTime>\n" % (time.time()) )
        xml_output.append("<ElapsedMinutes></ElapsedMinutes>\n")
        xml_output.append("</Testing>\n")
        xml_output.append("</Site>\n")

        with open(output_path, "w") as f:
            for line in xml_output:
                f.write(line)

        """
          /*
          <?xml version="1.0" encoding="UTF-8"?>
          <Site BuildName="Darwin-clang++"
                  BuildStamp="20121021-2300-Nightly"
                  Name="laphroaig.imp.fu-berlin.de"
                  Generator="ctest-2.8.9"
                  CompilerName="/usr/bin/clang++"
                  OSName="Mac OS X"
                  Hostname="laphroaig.imp.fu-berlin.de"
                  OSRelease="10.7.5"
                  OSVersion="11G63"
                  OSPlatform="x86_64"
                  Is64Bits="1"
                  VendorString="GenuineIntel"
                  VendorID="Intel Corporation"
                  FamilyID="6"
                  ModelID="37"
                  ProcessorCacheSize="32768"
                  NumberOfLogicalCPU="4"
                  NumberOfPhysicalCPU="2"
                  TotalVirtualMemory="512"
                  TotalPhysicalMemory="8192"
                  LogicalProcessorsPerPhysical="8"
                  ProcessorClockFrequency="2660"
          >
          <Testing>
          <StartDateTime>Oct 22 18:36 CEST</StartDateTime>
          <StartTestTime>1350923805</StartTestTime>
          <TestList>
            <Test>BinaryComposeFunctionAdapter_test</Test>
          </TestList>
          <Test Status="passed">
            <Name>BinaryComposeFunctionAdapter_test</Name>
            <Path>./source/TEST</Path>
            <FullName>./source/TEST/BinaryComposeFunctionAdapter_test</FullName>
            <FullCommandLine>/Users/aiche/dev/openms/openms-src/build/ninja/source/TEST/bin/BinaryComposeFunctionAdapter_test</FullCommandLine>
            <Results>
                    <NamedMeasurement type="numeric/double" name="Execution Time"><Value>0.469694</Value></NamedMeasurement>
                    <NamedMeasurement type="text/string" name="Completion Status"><Value>Completed</Value></NamedMeasurement>
                    <NamedMeasurement type="text/string" name="Command Line"><Value>/Users/aiche/dev/openms/openms-src/build/ninja/source/TEST/bin/BinaryComposeFunctionAdapter_test</Value></NamedMeasurement>
                    <Measurement>
                      <Value>
                      freier Text
                      </Value>
                    </Measurement>
            </Results>
          </Test>
              <EndDateTime>Oct 22 18:43 CEST</EndDateTime>
              <EndTestTime>1350924239</EndTestTime>
          <ElapsedMinutes>7.2</ElapsedMinutes></Testing>
        </Site>      
          */
                """


def checkPythonPxdHeader(src_path, bin_path, ignorefilename, pxds_out, print_pxd, output_format):
    """ Checks a set of doxygen xml file against a set of pxd header files

    For each C++ class found in the doxygen XML files, it tries to identify the
    corresponding pxd file. If a pxd file exists, it checks whether

    i)   all public functions, enums and attributes are wrapped in Python
    ii)  all void return types are correct in Python (these are not checked at
         compile time)
    iii) all fields of an enum are accessible from Python

    If it finds a method missing, the script suggests an addition and if a
    whole class is missing, the script writes suggestion .pxd file to a
    specified location (pxds_out).

    The output format can either be in text form (human readable) or in xml
    form which will try to overwrite the cdash Test.xml file to proivide an
    output to cdash. Please only specify xml output if in your binary
    you have executed "ctest -D Nightly" or similar.

    TODO also look at ./doc/doxygen/doxygen-error.log ?

    Make sure to build the doxygen xmls first with 
    $ make doc_xml

    """

    xml_output_path = os.path.join(bin_path, "doc", "xml_output")
    xml_files = glob.glob(xml_output_path + "/*.xml")
    print "Found %s doxygen xml files" % (len(xml_files))
    if len(xml_files) == 0:
        raise Exception("No doxygen files found in directory:\n%s,\n" % xml_output_path + \
                        "Please make sure you build the doxygen xmls (make dox_xml)\n" +\
                        "and that you specified the correct directory." )

    print "Creating pxd file map"
    pxd_file_matching = create_pxd_file_map(src_path)
    cnt = Counter()
    cnt.total = len(xml_files)
    ignorefile = IgnoreFile()
    if len(ignorefilename) > 0:
        ignorefile.load(ignorefilename)

    def pxd_text_printout(pxd_text, pxds_out, comp_name, print_pxd):
        if print_pxd: 
            print ""
            print pxd_text
        if len(pxds_out) > 0 and pxd_text is not None:
            fname =  os.path.join(pxds_out, "%s.pxd" % comp_name.split("::")[-1] )
            with open(fname, "w" ) as f:
                f.write(pxd_text)

    testresults = TestResultHandler()
    for class_cntr, f in enumerate(xml_files):
        dfile = DoxygenXMLFile(f)
        res = dfile.parse_doxygen()
        if dfile.parsing_error:
            # e.g. <computeroutput><bold></computeroutput>
            cnt.skipped_could_not_parse += 1
            msg = "Skip:: No-parse :: could not parse file %s with error %s" % (f, dfile.parsing_error_message)
            tres = TestResult(False, msg, name="%s_test" % f )
            testresults.append([ tres ])
            continue
        elif os.path.basename(f) == "index.xml":
            # Skip the index
            continue
        compound = res.get_compounddef()
        comp_name = compound.get_compoundname()
        if ignorefile.isNameIgnored(comp_name):
            msg = "Skip:: Ignored :: Class %s (file %s)" % (comp_name, f)
            tres = TestResult(True, msg, log_level=10, name="%s_test" % comp_name)
            testresults.append([ tres ])
            cnt.skipped_ignored += 1
            continue
        if compound.prot != "public":
            msg = "Skip:: Protected :: Compound %s is not public, skip" % (comp_name)
            tres = TestResult(True, msg, log_level=10, name="%s_test" % comp_name)
            testresults.append([ tres ])
            cnt.skipped_protected += 1
            continue
        file_location = dfile.getCompoundFileLocation()
        if file_location is None:
            msg = "Skip:: No-data :: there is no source file for %s" % f
            tres = TestResult(True, msg, log_level=10, name="%s_test" % comp_name)
            testresults.append([ tres ])
            cnt.skipped_no_location += 1
            continue
        openms_file = OpenMSSourceFile(file_location)
        maintainer = openms_file.getMaintainer()
        if dfile.isEmpty(True):
            msg = "Skip:: No-data :: File is empty (no section definitions found or only definitions found) in file %s" % f 
            tres = TestResult(True, msg, log_level=10, name="%s_test" % comp_name)
            tres.maintainer = maintainer
            testresults.append([ tres ])
            cnt.skipped_no_sections += 1
            continue
        if file_location in pxd_file_matching:
            pxdfiles = pxd_file_matching[file_location]
        else:
            msg = "Skip:: No-pxd :: No pxd file exists for Class %s (File %s) %s" % (comp_name, file_location, f)
            tres = TestResult(False, msg,  name="%s_test" % comp_name )
            tres.maintainer = maintainer
            testresults.append([ tres ])
            cnt.skipped_no_pxd_file += 1
            pxd_text = dfile.get_pxd_from_class(dfile, file_location, xml_output_path)
            pxd_text_printout(pxd_text, pxds_out, comp_name, print_pxd)
            continue
        try:
            pxd_class = PXDFile.parse_multiple_files(pxdfiles, comp_name)
            pxdfile = pxd_class.pxdfile
        except PXDFileParseError as e:
            # TODO specific exception
            msg = "Skip:: No-pxd :: " + e.message + "for %s (in pxd file %s)" % (comp_name, pxdfiles)
            tres = TestResult(False, msg,  name="%s_test" % comp_name )
            tres.maintainer = maintainer
            testresults.append([ tres ])
            cnt.skipped_no_pxd_match += 1
            pxd_text = dfile.get_pxd_from_class(dfile, file_location, xml_output_path)
            pxd_text_printout(pxd_text, pxds_out, comp_name, print_pxd)
            continue

        cnt.parsed += 1
        # Loop through all methods 
        classtestresults = []
        for method_cntr,mdef in enumerate(dfile.iterMemberDef()):

            if mdef.get_name() in ignorefile.getIgnoredMethods(comp_name):
                msg = "Ignore member function/attribute : %s %s %s " % (mdef.kind, mdef.prot, mdef.name)
                tres = TestResult(True, msg, log_level=10)
            else:
                tres = handle_member_definition(mdef, pxd_class, cnt)

            testname = "%s_%s_%s::%s" % (class_cntr, method_cntr, comp_name, mdef.name)
            testname = testname.replace("::", "_")
            testname = re.sub('[^a-zA-Z0-9_]+', '', testname)
            tres.comp_name = comp_name
            tres.file_location = file_location
            tres.pxdfile = pxdfile
            tres.maintainer = maintainer
            tres.name = testname
            classtestresults.append(tres)

        testresults.append(classtestresults)

    if output_format in ["text", "text-verbose", "text-quiet"]:
        for classtestresults in testresults:
            if len(classtestresults) > 1:
                t = classtestresults[0]
                lenfailed = len([t for t in classtestresults if not t.isPassed() ] )
                if lenfailed > 0:
                    print "== Test results for element %s - from Cpp file %s with maintainer %s and corresponding pxd file %s" % (
                        t.comp_name, t.file_location, t.maintainer, t.pxdfile)
            for tres in classtestresults:
                if not tres.isPassed():
                    print tres.message
                elif tres.log_level >= 10 and output_format in ["text", "text-verbose"]:
                    print tres.message
                elif tres.log_level >= 0 and output_format in ["text-verbose"]:
                    print tres.message

    elif output_format == "xml":

        # check if all files required to report in CDash are present
        tag_file = os.path.join(bin_path, "Testing", "TAG" ) 
        try: 
            # read the first line of tagfile (TAG) -> if it does not exist,
            # an IOError is thrown
            with open(tag_file) as f:
                ctestReportingPath = f.readline().strip()
                ctestReportingPath = os.path.join(bin_path, "Testing", ctestReportingPath)
                if not os.path.exists( ctestReportingPath ):
                    raise Exception("Missing directory at %s" % ( ctestReportingPath ) )
        except IOError:
            raise Exception("Missing nightly test information at %s" % (tag_file) )

        template_path = os.path.join(ctestReportingPath, "Test.xml" )
        testresults.to_cdash_xml(template_path, template_path)

    else:
        raise Exception("Unknown output format %s" % output_format)

    cnt.print_stats()
    cnt.print_skipping_reason()

def main(options):
    checkPythonPxdHeader(options.src_path, options.bin_path, options.ignorefile, options.pxds_out, options.print_pxd, options.output_format)

def handle_args():
    usage = "" 

    parser = argparse.ArgumentParser(description = usage )
    parser.add_argument("--bin_path", dest="bin_path", default=".", help="OpenMS build path")
    parser.add_argument("--src_path", dest="src_path", default=".", help="OpenMS source path")
    parser.add_argument("--ignore-file", dest="ignorefile", default="", help="Checker ignore file")
    parser.add_argument("--pxds-out", dest="pxds_out", default="", help="Folder to write pxd files")
    parser.add_argument("--output", dest="output_format", default="text", help="Output format (text or xml for ctest xml)")
    parser.add_argument('--print_pxd', action='store_true', default=False)
    #   print "Usage: checker.php <OpenMS src path> <OpenMS build path> [-u \"user name\"] [-t test] [options]\n";

    args = parser.parse_args(sys.argv[1:])
    return args

if __name__=="__main__":
    options = handle_args()
    main(options)


"""
offending doxygen lines that fail to parse:

include/source/APPLICATIONS/TOPP/IDRipper.C
  <B>NOTE: The meta value file origin is removed by the @p IDSplitter!!</B>
  generates
  the <computeroutput>IDSplitter!!</bold></computeroutput> 

doc/OpenMS_tutorial/OpenMS_Tutorial.doxygen
  \arg \c <b>[1]</b>:<A HREF="http://bieson.ub.uni-bielefeld.de/frontdoor.php?source_opus=1370"> 
  generates
  <listitem><para><computeroutput><bold></computeroutput>[1]</bold>:

"""

