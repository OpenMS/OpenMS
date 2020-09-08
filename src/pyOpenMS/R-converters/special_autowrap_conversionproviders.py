from __future__ import print_function
from autowrap.Code import Code
from autowrap.ConversionProvider import (TypeConverterBaseR,
                                         mangle,
                                         StdMapConverter)

class OpenMSDPosition2(TypeConverterBaseR):

    def get_base_types(self):
        return "DPosition2",

    def matches(self, cpp_type):
        return  not cpp_type.is_ptr

    def matching_python_type(self, cpp_type):
        return ""

    def type_check_expression(self, cpp_type, argument_var):
        return "is.vector(%s) && length(%s) == 2 && (is_scalar_integer(%s[[1]]) || is_scalar_double(%s[[1]])) && (is_scalar_integer(%s[[2]]) || is_scalar_double(%s[[2]]))" % (argument_var, argument_var, argument_var,argument_var, argument_var, argument_var)

    def input_conversion(self, cpp_type, argument_var, arg_num):
        cr_ref = False
        dp = "dp_%s" % arg_num
        code = Code().add("""
            |$dp <- r_to_py($argument_var)
        """, locals())
        cleanup = ""
        if cpp_type.is_ref:
            cr_ref = True
            cleanup = Code().add("""
            |by_ref${arg_num} = py_to_r($dp)
            """, locals())
        call_as = dp
        return code, call_as, cleanup, cr_ref

    def output_conversion(self, cpp_type, input_r_var, output_r_var):
        # this one is slow as it uses construction of python type DataValue for
        # delegating conversion to this type, which reduces code below:
        return Code().add("""
                    |$output_r_var = $input_r_var
                """, locals())


class OpenMSDPosition2Vector(TypeConverterBaseR):

    def get_base_types(self):
        return "libcpp_vector",

    def matches(self, cpp_type):
        inner_t, = cpp_type.template_args
        return inner_t == "DPosition2"

    def matching_python_type(self, cpp_type):
        return "np.ndarray[np.float32_t,ndim=2]"

    def type_check_expression(self, cpp_type, argument_var):
        # Reticulate converts matrix/array to numpy array.
        # for this type matrix is better suited.
        return "is.matrix(%s) && NROW(%s) == 2 && is_double(%s[1,]) && is_double(%s[2,])" % (argument_var,argument_var,argument_var,argument_var)

    def input_conversion(self, cpp_type, argument_var, arg_num):
        dp = "dp_%s" % arg_num
        cr_ref = False
        code = Code().add("""
            |$dp <- r_to_py($argument_var)
        """, locals())
        cleanup = ""
        if cpp_type.is_ref:
            cr_ref = True
            cleanup = Code().add("""
            |byref_${arg_num} <- py_to_r($dp)
            """, locals())
        call_as = dp
        return code, call_as, cleanup, cr_ref

    def output_conversion(self, cpp_type, input_r_var, output_r_var):
        return Code().add("""
         |$output_r_var <- $input_r_var
         """, locals())


class OpenMSDataValue(TypeConverterBaseR):

    def get_base_types(self):
        return "DataValue",

    def matches(self, cpp_type):
        return  not cpp_type.is_ptr and not cpp_type.is_ref

    def matching_python_type(self, cpp_type):
        return ""

    def type_check_expression(self, cpp_type, argument_var):
        return "is_scalar_integer(%s) || is_scalar_double(%s) || is.vector(%s) || is_scalar_character(%s)" % (argument_var,argument_var,argument_var,argument_var)

    def input_conversion(self, cpp_type, argument_var, arg_num):
        cr_ref = False
        arg_conv = self.converters.get(cpp_type)
        call_as = "%s" % argument_var
        return "", call_as, "", cr_ref

    def output_conversion(self, cpp_type, input_cpp_var, output_py_var):
        # this one is slow as it uses construction of python type DataValue for
        # delegating conversion to this type, which reduces code below:
        return Code().add("""
                    |if(length($input_cpp_var) > 1){
                    |   if( all(sapply($input_cpp_var,function(x) class(x)[1]=="python.builtin.bytes")) ){
                    |   $output_py_var = sapply($input_cpp_var,as.character)
                    |   }
                    |} else if (length($input_cpp_var) == 1){
                    |   $output_py_var = $input_cpp_var
                    |} else if (is.null($input_cpp_var)){
                    |   $output_py_var = NULL
                    |}  
                """, locals())


class OpenMSStringConverter(TypeConverterBaseR):

    def get_base_types(self):
        return "String",

    def matches(self, cpp_type):
        return  not cpp_type.is_ptr

    def matching_python_type(self, cpp_type):
        # We allow bytes, unicode, str and String
        return ""

    def type_check_expression(self, cpp_type, argument_var):
        # Need to treat ptr and reference differently as these may be modified
        # and the results needs to be available in Python
        if (cpp_type.is_ptr or cpp_type.is_ref) and not cpp_type.is_const:
            pass

        # Allow conversion from unicode str, bytes and OpenMS::String
        return "(is.R6(%s) && class(%s)[1]==\"String\") || is_scalar_character(%s)" % (argument_var, argument_var, argument_var)

    def input_conversion(self, cpp_type, argument_var, arg_num):

        # See ./src/pyOpenMS/addons/ADD_TO_FIRST.pyx for declaration of convString
        cr_ref = False
        call_as = "r_to_py(%s)" % argument_var
        cleanup = ""
        code = ""
        return code, call_as, cleanup, cr_ref

    def output_conversion(self, cpp_type, input_cpp_var, output_py_var):

        # See ./src/pyOpenMS/addons/ADD_TO_FIRST.pyx for declaration of convOutputString
        return "%s = %s" % (output_py_var, input_cpp_var)

class AbstractOpenMSListConverter(TypeConverterBaseR):

    def __init__(self):
        # mark as abstract:
        raise NotImplementedError()

    def get_base_types(self):
        return self.openms_type,

    def matches(self, cpp_type):
        return not cpp_type.is_ptr

    def matching_python_type(self, cpp_type):
        return "list"

    def type_check_expression(self, cpp_type, argument_var):
        if self.inner_py_type == "int":
            check_inner_type = "is_scalar_integer"
        elif self.inner_py_type == "float":
            check_inner_type = "is_scalar_double"
        return "(is_list(%s) || is_vector(%s) && (length(%s) > 1 || length(%s)==0)) && all(sapply(%s, %s))" % (argument_var, argument_var, argument_var, argument_var, argument_var,check_inner_type)

    def input_conversion(self, cpp_type, argument_var, arg_num):
        cr_ref = False
        temp_var = "v%d" % arg_num
        t = self.inner_cpp_type
        ltype = self.openms_type
        code = Code().add("""
                |$temp_var = r_to_py($argument_var)
                |type_check = ifelse(is_list($argument_var),"list","vector")
                """, locals())
        cleanup = ""
        if cpp_type.is_ref:
            cr_ref = True
            cleanup_code = Code().add("""
                    |if (type_check == "list"){
                    |   byref_${arg_num} <- as.list(py_to_r($temp_var))
                    |}
                    |else{
                    |   byref_${arg_num} <- py_to_r($temp_var)
                    |}
                    """, locals())
        call_as = "%s" % temp_var
        return code, call_as, cleanup, cr_ref

    def output_conversion(self, cpp_type, input_cpp_var, output_py_var):
        t = self.inner_cpp_type
        code = Code().add("""
            |if(length($input_cpp_var)==1){
            |   $output_py_var = as.list($input_cpp_var)
            |} else { $output_py_var = $input_cpp_var }
            """, locals())
        return code

class OpenMSIntListConverter(AbstractOpenMSListConverter):

    openms_type = "IntList"
    inner_py_type = "int"
    inner_cpp_type = "int"
    # mark as non abstract:
    def __init__(self):
        pass

class OpenMSDoubleListConverter(AbstractOpenMSListConverter):

    openms_type = "DoubleList"
    inner_py_type = "float"
    inner_cpp_type = "double"
    # mark as non abstract:
    def __init__(self):
        pass

class StdVectorStringConverter(TypeConverterBaseR):

    def get_base_types(self):
        return "libcpp_vector",

    def matches(self, cpp_type):
        inner_t, = cpp_type.template_args
        return inner_t == "String"

    def matching_python_type(self, cpp_type):
        return "list"

    def type_check_expression(self, cpp_type, arg_var):
        return Code().add("""
          |is.vector($arg_var) && all(sapply($arg_var,is_scalar_character))
          """, locals()).render()

    def input_conversion(self, cpp_type, argument_var, arg_num):
        cr_ref = False
        temp_var = "v%d" % arg_num
        code = Code().add("""
            |$temp_var = r_to_py(lapply($argument_var,function(a) py_builtin$$bytes(a,'utf-8')))
            """, locals())
        if cpp_type.is_ref:
            cr_ref = True
            cleanup_code = Code().add("""
                |byref_${arg_num} <- sapply(py_to_r($temp_var),as.character)
                """, locals())
        else:
            cleanup_code = ""
        return code, "%s" % temp_var, cleanup_code, cr_ref

    def call_method(self, res_type, cy_call_str):
        return "py_ans = %s" % (cy_call_str)

    def output_conversion(self, cpp_type, input_cpp_var, output_py_var):

        if cpp_type.is_ptr:
            raise AssertionError()

        code = Code().add("""
            |$output_py_var = sapply($input_cpp_var,as.character)
            """, locals())
        return code

class OpenMSStringListConverter(StdVectorStringConverter):
    """
    StringList is now a std::vector<String> so we can directly
    inherit all methods from StdVectorStringConverter but need 
    to add the matching rules for StringList.
    """
    openms_type = "StringList"
    inner_py_type = "bytes"
    inner_cpp_type = "libcpp_string"
    # mark as non abstract:
    def __init__(self):
        pass

    def get_base_types(self):
        return self.openms_type,

    def matches(self, cpp_type):
        return not cpp_type.is_ptr

    def matching_python_type(self, cpp_type):
        return "list"
    
    def type_check_expression(self, cpp_type, arg_var):
        return Code().add("""
        |(is_vector($arg_var) && (length($arg_var) > 1 || length($arg_var)==0) || is_list($arg_var)) && all(sapply($arg_var,is_scalar_character))
        """, locals()).render()

    def output_conversion(self, cpp_type, input_cpp_var, output_py_var):
    
        if cpp_type.is_ptr:
            raise AssertionError()

        code = Code().add("""
            |   if(length($input_cpp_var) == 1){
            |       $output_py_var <- lapply($input_cpp_var, as.character)
            |   } else { $output_py_var <- sapply($input_cpp_var,as.character) } 
            """, locals())
        return code
    
class StdSetStringConverter(TypeConverterBaseR):

    def get_base_types(self):
        return "libcpp_set",

    def matches(self, cpp_type):
        inner_t, = cpp_type.template_args
        return inner_t == "String"

    def matching_python_type(self, cpp_type):
        return "set"

    def type_check_expression(self, cpp_type, arg_var):
        return Code().add("""
          |is.vector($arg_var) && all(sapply($arg_var,is_scalar_character)) && !any(duplicated($arg_var) == T)
          """, locals()).render()

    def input_conversion(self, cpp_type, argument_var, arg_num):
        cr_ref = False
        temp_var = "v%d" % arg_num
        code = Code().add("""
            |$temp_var = py_builtin$$set(modify_depth($argument_var,1,function(a) py_builtin$$bytes(a,'utf-8')))
            """, locals())
        if cpp_type.is_ref:
            cr_ref = True
            cleanup_code = Code().add("""
                |byref_${arg_num} <- modify_depth(py_to_r(py_builtin$$list($temp_var)),1,as.character)
                """, locals())
        else:
            cleanup_code = ""
        return code, "%s" % temp_var, cleanup_code, cr_ref

    def call_method(self, res_type, cy_call_str):
        return "py_ans = %s" % (cy_call_str)


    def output_conversion(self, cpp_type, input_cpp_var, output_py_var):

        if cpp_type.is_ptr:
            raise AssertionError()

        code = Code().add("""
            |$output_py_var = modify_depth(py_to_r(py_builtin$$list($input_cpp_var)),1,as.character)
            """, locals())
        return code

class OpenMSMapConverter(StdMapConverter):

    def get_base_types(self):
        return "Map",

    def matches(self, cpp_type):
        return True

    def matching_python_type(self, cpp_type):
        return "dict"

    def type_check_expression(self, cpp_type, arg_var):
        tt_key, tt_value = cpp_type.template_args
        inner_conv_1 = self.converters.get(tt_key)
        inner_conv_2 = self.converters.get(tt_value)

        if inner_conv_1 is None:
            raise Exception("arg type %s not supported" % tt_key)
        if inner_conv_2 is None:
            raise Exception("arg type %s not supported" % tt_value)

        inner_check_1 = inner_conv_1.type_check_expression(tt_key, "k")
        inner_check_2 = inner_conv_2.type_check_expression(tt_value, "v")

        return """
          is.environment(%s) && identical(parent.env(%s), asNamespace("collections")) && identical(strsplit(capture.output(%s$print())," ")[[1]][1], "dict")
          && all(sapply(%s$keys(),function(k) %s))
          && all(sapply(%s$values(),function(v) %s))
          """ % (arg_var,arg_var,arg_var,arg_var,inner_check_1,arg_var,inner_check_2)

    def input_conversion(self, cpp_type, argument_var, arg_num):
        cr_ref = False
        tt_key, tt_value = cpp_type.template_args
        temp_var = "v%d" % arg_num

        cy_tt_key = self.converters.cython_type(tt_key)
        cy_tt_value = self.converters.cython_type(tt_value)

        if cy_tt_key.is_enum:
            key_conv = "%s$keys()" % argument_var
        elif tt_key.base_type in self.converters.names_to_wrap:
            raise Exception("can not handle wrapped classes as keys in map")
        else:
            key_conv = "%s$keys()" % argument_var

        if cy_tt_value.is_enum:
            value_conv = "%s$values()" % argument_var
        elif tt_value.base_type in self.converters.names_to_wrap:
            value_conv = "lapply(%s$values(),function(v) r_to_py(v))" % argument_var
        else:
            value_conv = "%s$values()" % argument_var

        code = Code().add("""
            |$temp_var = py_dict(${key_conv},${value_conv})
            """, locals())

        if cpp_type.is_ref:
            cr_ref = True
            if cy_tt_key.is_enum:
                key_conv = "py_to_r(py_builtin$list(%s$keys()))" % (temp_var)
            elif tt_key.base_type in self.converters.names_to_wrap:
                raise Exception("can not handle wrapped classes as keys in map")
            else:
                key_conv = "py_to_r(py_builtin$list(%s$keys()))" % (temp_var)

            if not cy_tt_value.is_enum and tt_value.base_type in self.converters.names_to_wrap:
                cy_tt = tt_value.base_type
                value_conv = "py_to_r(py_builtin$list(%s$values()))" % (temp_var)
                cleanup_code = Code().add("""
                    |byref_${arg_num} <- collections::dict(lapply(${value_conv},function(v) ${cy_tt}$$new(v)), ${key_conv})
                    """, locals())
            else:
                value_conv = "py_to_r(py_builtin$list(%s$values()))" % (temp_var)
                cleanup_code = Code().add("""
                    |byref_${arg_num} <- collections::dict(${value_conv},${key_conv})
                    """, locals())
        else:
            cleanup_code = ""

        return code, "%s" % temp_var, cleanup_code, cr_ref

    def call_method(self, res_type, cy_call_str):
        return "py_ans = %s" % (cy_call_str)

    def output_conversion(self, cpp_type, input_cpp_var, output_py_var):

        if cpp_type.is_ptr:
            raise AssertionError()

        tt_key, tt_value = cpp_type.template_args
        cy_tt_key = self.converters.cython_type(tt_key)
        cy_tt_value = self.converters.cython_type(tt_value)

        if not cy_tt_key.is_enum and tt_key.base_type in self.converters.names_to_wrap:
            raise Exception("can not handle wrapped classes as keys in map")
        else:
            key_conv = "py_to_r(py_builtin$list(%s$keys()))" % (input_cpp_var)

        if not cy_tt_value.is_enum and tt_value.base_type in self.converters.names_to_wrap:
            cy_tt = tt_value.base_type
            value_conv = "py_to_r(py_builtin$list(%s$values()))" % (input_cpp_var)
            code = Code().add("""
                |$output_py_var = collections::dict(lapply(${value_conv},function(v) ${cy_tt}$$new(v)), ${key_conv})
                """, locals())
            return code
        else:
            value_conv = "py_to_r(py_builtin$list(%s$values()))" % (input_cpp_var)
            code = Code().add("""
                |$output_py_var = collections::dict(${value_conv}, ${key_conv})
                """, locals())
            return code


import time

class CVTermMapConverter(TypeConverterBaseR):

    def get_base_types(self):
        return "Map",

    def matches(self, cpp_type):
        return str(cpp_type) == "Map[String,libcpp_vector[CVTerm]]" \
           or  str(cpp_type) == "Map[String,libcpp_vector[CVTerm]] &"

    def matching_python_type(self, cpp_type):
        return "dict"

    def type_check_expression(self, cpp_type, arg_var):
        return """
          is.environment(%s) && identical(parent.env(%s), asNamespace("collections")) && identical(strsplit(capture.output(%s$print())," ")[[1]][1], "dict")
          && all(sapply(%s$keys(),is_scalar_character))
          && all(sapply(%s$values(), function(v) is.vector(v) && sapply(v, function(v1) is.R6(v1) && class(v1)[1] == "CVTerm")))
          """ % (arg_var,arg_var,arg_var,arg_var,arg_var)

    def input_conversion(self, cpp_type, argument_var, arg_num):
        cr_ref = False
        map_name = "map_%d" % arg_num

        key_conv = "modify_depth(%s$keys(),1,function(i) py_builtin$bytes(i,'utf-8'))" % (argument_var)
        value_conv = "modify_depth(%s$values(),2,function(v) r_to_py(v))" % (argument_var)
        code = Code().add("""
                |$map_name <- py_dict(${value_conv},${key_conv})
                """, locals())

        if cpp_type.is_ref:
            cr_ref = True
            key_conv = "py_to_r(py_builtin$list(%s$keys()))" % (map_name)
            value_conv = "py_to_r(py_builtin$list(%s$values()))" % (map_name)
            cleanup_code = Code().add("""
                |byref_${arg_num} <- collections::dict(lapply(${value_conv},function(v) CVTerm$$new(v)), lapply(${key_conv},as.character))
                """, locals())
        else:
            cleanup_code = ""
        return code, map_name, cleanup_code, cr_ref


    def call_method(self, res_type, cy_call_str):
        return "py_ans = %s" % (cy_call_str)

    def output_conversion(self, cpp_type, input_cpp_var, output_py_var):

        key_conv = "py_to_r(py_builtin$list(%s$keys()))" % (input_cpp_var)
        value_conv = "py_to_r(py_builtin$list(%s$values()))" % (input_cpp_var)
        code = Code().add("""
            |$output_py_var <- collections::dict(lapply(${value_conv},function(v) CVTerm$$new(v)), lapply(${key_conv},as.character))
            """, locals())

        return code


