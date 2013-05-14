from autowrap.Code import Code
from autowrap.ConversionProvider import TypeConverterBase, mangle


class OpenMSStringConverter(TypeConverterBase):

    def get_base_types(self):
        return "String",

    def matches(self, cpp_type):
        return  not cpp_type.is_ptr and not cpp_type.is_ref

    def matching_python_type(self, cpp_type):
        return "bytes"

    def type_check_expression(self, cpp_type, argument_var):
        return "isinstance(%s, bytes)" % (argument_var,)

    def input_conversion(self, cpp_type, argument_var, arg_num):
        # here we inject special behavoir for testing if this converter
        # was called !
        call_as = "(_String(<char *>%s))" % argument_var
        code = cleanup = ""
        return code, call_as, cleanup

    def output_conversion(self, cpp_type, input_cpp_var, output_py_var):
        return "%s = _cast_const_away(<char*>%s.c_str())" % (output_py_var,
                                                             input_cpp_var)


class AbstractOpenMSListConverter(TypeConverterBase):

    def __init__(self):
        # mark as abstract:
        raise NotImplementedError()

    def get_base_types(self):
        return self.openms_type,

    def matches(self, cpp_type):
        return  not cpp_type.is_ptr and not cpp_type.is_ref

    def matching_python_type(self, cpp_type):
        return "list"

    def type_check_expression(self, cpp_type, argument_var):
        return\
            "isinstance(%s, list) and all(isinstance(li, %s) for li in %s)"\
            % (argument_var, self.inner_py_type, argument_var)

    def input_conversion(self, cpp_type, argument_var, arg_num):
        temp_var = "v%d" % arg_num
        t = self.inner_cpp_type
        ltype = self.openms_type
        code = Code().add("""
                |cdef libcpp_vector[$t] _$temp_var = $argument_var
                |cdef _$ltype $temp_var = _$ltype(_$temp_var)
                """, locals())
        cleanup = ""
        # here we inject special behavoir for testing if this converter
        # was called !
        call_as = "(%s)" % temp_var
        return code, call_as, cleanup

    def output_conversion(self, cpp_type, input_cpp_var, output_py_var):
        t = self.inner_cpp_type
        code = Code().add("""
            |$output_py_var = []
            |cdef int i, n
            |n = $input_cpp_var.size()
            |for i in range(n):
            |    $output_py_var.append(<$t>$input_cpp_var.at(i))
            """, locals())
        return code


class OpenMSIntListConverter(AbstractOpenMSListConverter):

    openms_type = "IntList"
    inner_py_type = "int"
    inner_cpp_type = "int"
    # mark as non abstract:
    def __init__(self):
        pass

class OpenMSStringListConverter(AbstractOpenMSListConverter):

    openms_type = "StringList"
    inner_py_type = "bytes"
    inner_cpp_type = "libcpp_string"
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

class StdVectorStringConverter(TypeConverterBase):

    def get_base_types(self):
        return "libcpp_vector",

    def matches(self, cpp_type):
        inner_t, = cpp_type.template_args
        return inner_t == "String"

    def matching_python_type(self, cpp_type):
        return "list"

    def type_check_expression(self, cpp_type, arg_var):
        return Code().add("""
          |isinstance($arg_var, list) and all(isinstance(i, bytes) for i in
          + $arg_var)
          """, locals()).render()

    def input_conversion(self, cpp_type, argument_var, arg_num):
        temp_var = "v%d" % arg_num
        item = "item%d" % arg_num
        code = Code().add("""
            |cdef libcpp_vector[_String] * $temp_var = new libcpp_vector[_String]()
            |cdef bytes $item
            |for $item in $argument_var:
            |   $temp_var.push_back(_String(<char *>$item))
            """, locals())
        if cpp_type.is_ref:
            cleanup_code = Code().add("""
                |cdef replace = []
                |cdef libcpp_vector[_String].iterator it = $temp_var.begin()
                |while it != $temp_var.end():
                |   replace.append(<char*>deref(it).c_str())
                |   inc(it)
                |$argument_var[:] = replace
                |del $temp_var
                """, locals())
        else:
            cleanup_code = "del %s" % temp_var
        return code, "deref(%s)" % temp_var, cleanup_code

    def call_method(self, res_type, cy_call_str):
        return "_r = %s" % (cy_call_str)


    def output_conversion(self, cpp_type, input_cpp_var, output_py_var):

        assert not cpp_type.is_ptr
        it = mangle("it_" + input_cpp_var)
        item = mangle("item_" + output_py_var)
        code = Code().add("""
            |$output_py_var = []
            |cdef libcpp_vector[_String].iterator $it = $input_cpp_var.begin()
            |while $it != $input_cpp_var.end():
            |   $output_py_var.append(<char*>deref($it).c_str())
            |   inc($it)
            """, locals())
        return code

class StdSetStringConverter(TypeConverterBase):

    def get_base_types(self):
        return "libcpp_set",

    def matches(self, cpp_type):
        inner_t, = cpp_type.template_args
        return inner_t == "String"

    def matching_python_type(self, cpp_type):
        return "set"

    def type_check_expression(self, cpp_type, arg_var):
        return Code().add("""
          |isinstance($arg_var, set) and all(isinstance(i, bytes) for i in
          + $arg_var)
          """, locals()).render()

    def input_conversion(self, cpp_type, argument_var, arg_num):
        temp_var = "v%d" % arg_num
        item = "item%d" % arg_num
        code = Code().add("""
            |cdef libcpp_set[_String] * $temp_var = new libcpp_set[_String]()
            |cdef bytes $item
            |for $item in $argument_var:
            |   $temp_var.insert(_String(<char *>$item))
            """, locals())
        if cpp_type.is_ref:
            cleanup_code = Code().add("""
                |cdef replace = set()
                |cdef libcpp_set[_String].iterator it = $temp_var.begin()
                |while it != $temp_var.end():
                |   replace.add(<char*>deref(it).c_str())
                |   inc(it)
                |$argument_var.clear()
                |$argument_var.update(replace)
                |del $temp_var
                """, locals())
        else:
            cleanup_code = "del %s" % temp_var
        return code, "deref(%s)" % temp_var, cleanup_code

    def call_method(self, res_type, cy_call_str):
        return "_r = %s" % (cy_call_str)


    def output_conversion(self, cpp_type, input_cpp_var, output_py_var):

        assert not cpp_type.is_ptr
        it = mangle("it_" + input_cpp_var)
        item = mangle("item_" + output_py_var)
        code = Code().add("""
            |$output_py_var = set()
            |cdef libcpp_set[_String].iterator $it = $input_cpp_var.begin()
            |while $it != $input_cpp_var.end():
            |   $output_py_var.add(<char*>deref($it).c_str())
            |   inc($it)
            """, locals())
        return code




def register_all():
    from autowrap.ConversionProvider import  special_converters
    special_converters.append(OpenMSStringConverter())
    special_converters.append(OpenMSStringListConverter())
    special_converters.append(OpenMSIntListConverter())
    special_converters.append(OpenMSDoubleListConverter())
    special_converters.append(StdVectorStringConverter())
    special_converters.append(StdSetStringConverter())

