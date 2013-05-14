from autowrap.Code import Code
from autowrap.ConversionProvider import (TypeConverterBase,
                                         mangle,
                                         StdMapConverter)


class OpenMSDPosition2(TypeConverterBase):

    def get_base_types(self):
        return "DPosition2",

    def matches(self, cpp_type):
        return  not cpp_type.is_ptr

    def matching_python_type(self, cpp_type):
        return ""

    def type_check_expression(self, cpp_type, argument_var):
        return "len(%s) == 2 and isinstance(%s[0], (int, float)) "\
                "and isinstance(%s[1], (int, float))"\
                % (argument_var, argument_var, argument_var)

    def input_conversion(self, cpp_type, argument_var, arg_num):
        dp = "_dp_%s" % arg_num
        code = Code().add("""
            |cdef _DPosition2 $dp
            |$dp[0] = <float>$argument_var[0]
            |$dp[1] = <float>$argument_var[1]
        """, locals())
        cleanup = ""
        if cpp_type.is_ref:
            cleanup = Code().add("""
            |$cpp_type[0] = $dp[0]
            |$cpp_type[1] = $dp[1]
            """, locals())
        call_as = dp
        return code, call_as, cleanup

    def output_conversion(self, cpp_type, input_cpp_var, output_py_var):
        # this one is slow as it uses construction of python type DataValue for
        # delegating conversion to this type, which reduces code below:
        return Code().add("""
                    |$output_py_var = [$input_cpp_var[0], $input_cpp_var[1]]
                """, locals())


class OpenMSDPosition2Vector(TypeConverterBase):

    def get_base_types(self):
        return "libcpp_vector",

    def matches(self, cpp_type):
        inner_t, = cpp_type.template_args
        return inner_t == "DPosition2"

    def matching_python_type(self, cpp_type):
        return "np.ndarray[np.float32_t,ndim=2]"

    def type_check_expression(self, cpp_type, argument_var):
        return "%s.shape[1] == 2" % argument_var

    def input_conversion(self, cpp_type, argument_var, arg_num):
        dp = "_dp_%s" % arg_num
        vec ="_dp_vec_%s" % arg_num
        ii = "_dp_ii_%s" % arg_num
        N = "_dp_N_%s" % arg_num
        code = Code().add("""
            |cdef libcpp_vector[_DPosition2] $vec
            |cdef _DPosition2 $dp
            |cdef int $ii
            |cdef int $N = $argument_var.shape[0]
            |for $ii in range($N):
            |    $dp[0] = $argument_var[$ii,0]
            |    $dp[1] = $argument_var[$ii,1]
            |    $vec.push_back($dp)
        """, locals())
        cleanup = ""
        if cpp_type.is_ref:
            it = "_dp_it_%s" % arg_num
            n = "_dp_n_%s" % arg_num
            cleanup = Code().add("""
            |$n = $vec.size()
            |$argument_var.resize(($n,2))

            |$n = $vec.size()
            |cdef libcpp_vector[_DPosition2].iterator $it = $vec.begin()
            |$ii = 0
            |while $it != $vec.end():
            |     $argument_var[$ii, 0] = deref($it)[0]
            |     $argument_var[$ii, 1] = deref($it)[1]
            |     inc($it)
            |     $ii += 1

            """, locals())
        call_as = vec
        return code, call_as, cleanup

    def output_conversion(self, cpp_type, input_cpp_var, output_py_var):
        # this one is slow as it uses construction of python type DataValue for
        # delegating conversion to this type, which reduces code below:
        it = "_out_it_dpos_vec"
        n  = "_out_n_dpos_vec"
        ii = "_out_ii_dpos_vec"
        return Code().add("""
         |cdef int $n = $input_cpp_var.size()
         |cdef $output_py_var = np.zeros([$n,2], dtype=np.float32)
         |cdef libcpp_vector[_DPosition2].iterator $it = $input_cpp_var.begin()
         |cdef int $ii = 0
         |while $it != $input_cpp_var.end():
         |     $output_py_var[$ii, 0] = deref($it)[0]
         |     $output_py_var[$ii, 1] = deref($it)[1]
         |     inc($it)
         |     $ii += 1
         """, locals())


class OpenMSDataValue(TypeConverterBase):

    def get_base_types(self):
        return "DataValue",

    def matches(self, cpp_type):
        return  not cpp_type.is_ptr and not cpp_type.is_ref

    def matching_python_type(self, cpp_type):
        return ""

    def type_check_expression(self, cpp_type, argument_var):
        return "isinstance(%s, (int, long, float, list, str))" % argument_var

    def input_conversion(self, cpp_type, argument_var, arg_num):
        call_as = "deref(DataValue(%s).inst.get())" % argument_var
        return "", call_as, ""

    def output_conversion(self, cpp_type, input_cpp_var, output_py_var):
        # this one is slow as it uses construction of python type DataValue for
        # delegating conversion to this type, which reduces code below:
        return Code().add("""
                    |cdef DataValue _vvv = DataValue.__new__(DataValue)
                    |_vvv.inst = shared_ptr[_DataValue](new _DataValue($input_cpp_var))
                    |cdef int _value_type_x = $input_cpp_var.valueType()
                    |if _value_type_x == DataType.STRING_VALUE:
                    |    $output_py_var = _vvv.toString()
                    |elif _value_type_x == DataType.INT_VALUE:
                    |    $output_py_var = _vvv.toInt()
                    |elif _value_type_x == DataType.DOUBLE_VALUE:
                    |    $output_py_var = _vvv.toDouble()
                    |elif _value_type_x == DataType.INT_LIST:
                    |    $output_py_var = _vvv.toIntList()
                    |elif _value_type_x == DataType.DOUBLE_LIST:
                    |    $output_py_var = _vvv.toDoubleList()
                    |elif _value_type_x == DataType.STRING_LIST:
                    |    $output_py_var = _vvv.toStringList()
                """, locals())


class OpenMSStringConverter(TypeConverterBase):

    def get_base_types(self):
        return "String",

    def matches(self, cpp_type):
        return  not cpp_type.is_ptr

    def matching_python_type(self, cpp_type):
        return "String" if cpp_type.is_ref else "bytes"

    def type_check_expression(self, cpp_type, argument_var):
        if cpp_type.is_ref:
            return "isinstance(%s, String)" % (argument_var,)
        return "isinstance(%s, bytes)" % (argument_var,)

    def input_conversion(self, cpp_type, argument_var, arg_num):
        # here we inject special behavoir for testing if this converter
        # was called !
        if cpp_type.is_ref:
            call_as = "deref(%s.inst.get())" % argument_var
        else:
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
        return  not cpp_type.is_ptr

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
        if cpp_type.is_ref:
            cleanup_code = Code().add("""
                    |replace = []
                    |cdef int i, n
                    |n = $temp_var.size()
                    |for i in range(n):
                    |    replace.append(<$t>$temp_var.at(i))
                    |$argument_var[:] = replace
                    """, locals())
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
        assert inner_conv_1 is not None, "arg type %s not supported" % tt_key
        assert inner_conv_2 is not None, "arg type %s not supported" % tt_value

        inner_check_1 = inner_conv_1.type_check_expression(tt_key, "k")
        inner_check_2 = inner_conv_2.type_check_expression(tt_value, "v")

        return Code().add("""
          |isinstance($arg_var, dict)
          + and all($inner_check_1 for k in $arg_var.keys())
          + and all($inner_check_2 for v in $arg_var.values())
          """, locals()).render()

    def input_conversion(self, cpp_type, argument_var, arg_num):
        tt_key, tt_value = cpp_type.template_args
        temp_var = "v%d" % arg_num

        cy_tt_key = self.converters.cython_type(tt_key)
        cy_tt_value = self.converters.cython_type(tt_value)

        if cy_tt_key.is_enum:
            key_conv = "<%s> key" % cy_tt_key
        elif tt_key.base_type in self.converters.names_to_wrap:
            raise Exception("can not handle wrapped classes as keys in map")
        else:
            key_conv = "<%s> key" % cy_tt_key

        if cy_tt_value.is_enum:
            value_conv = "<%s> value" % cy_tt_value
        elif tt_value.base_type in self.converters.names_to_wrap:
            value_conv = "deref((<%s>value).inst.get())" % tt_value.base_type
        else:
            value_conv = "<%s> value" % cy_tt_value

        code = Code().add("""
            |cdef _Map[$cy_tt_key, $cy_tt_value] * $temp_var = new
            + _Map[$cy_tt_key, $cy_tt_value]()

            |for key, value in $argument_var.items():
            |   deref($temp_var)[$key_conv] = $value_conv
            """, locals())

        if cpp_type.is_ref:
            replace = mangle("replace_" + argument_var)
            it = mangle("it_" + argument_var)

            if cy_tt_key.is_enum:
                key_conv = "<%s> deref(%s).first" % (cy_tt_key, it)
            elif tt_key.base_type in self.converters.names_to_wrap:
                raise Exception("can not handle wrapped classes as keys in map")
            else:
                key_conv = "<%s> deref(%s).first" % (cy_tt_key, it)

            if not cy_tt_value.is_enum and tt_value.base_type in self.converters.names_to_wrap:
                cy_tt = tt_value.base_type
                item = mangle("item_" + argument_var)
                cleanup_code = Code().add("""
                    |cdef $replace = dict()
                    |cdef _Map[$cy_tt_key, $cy_tt_value].iterator $it = $temp_var.begin()
                    |cdef $cy_tt $item
                    |while $it != $temp_var.end():
                    |   $item = $cy_tt.__new__($cy_tt)
                    |   $item.inst = shared_ptr[$cy_tt_value](new $cy_tt_value((deref($it)).second))
                    |   $replace[$key_conv] = $item
                    |   inc($it)
                    |$argument_var.clear()
                    |$argument_var.update($replace)
                    |del $temp_var
                    """, locals())
            else:
                value_conv = "<%s> deref(%s).second" % (cy_tt_value, it)
                cleanup_code = Code().add("""
                    |cdef $replace = dict()
                    |cdef _Map[$cy_tt_key, $cy_tt_value].iterator $it = $temp_var.begin()
                    |while $it != $temp_var.end():
                    |   $replace[$key_conv] = $value_conv
                    |   inc($it)
                    |$argument_var.clear()
                    |$argument_var.update($replace)
                    |del $temp_var
                    """, locals())
        else:
            cleanup_code = "del %s" % temp_var

        return code, "deref(%s)" % temp_var, cleanup_code

    def call_method(self, res_type, cy_call_str):
        return "_r = %s" % (cy_call_str)

    def output_conversion(self, cpp_type, input_cpp_var, output_py_var):

        assert not cpp_type.is_ptr

        tt_key, tt_value = cpp_type.template_args
        cy_tt_key = self.converters.cython_type(tt_key)
        cy_tt_value = self.converters.cython_type(tt_value)

        it = mangle("it_" + input_cpp_var)

        if not cy_tt_key.is_enum and tt_key.base_type in self.converters.names_to_wrap:
            raise Exception("can not handle wrapped classes as keys in map")
        else:
            key_conv = "<%s>(deref(%s).first)" % (cy_tt_key, it)

        if not cy_tt_value.is_enum and tt_value.base_type in self.converters.names_to_wrap:
            cy_tt = tt_value.base_type
            item = mangle("item_" + output_py_var)
            code = Code().add("""
                |$output_py_var = dict()
                |cdef _Map[$cy_tt_key, $cy_tt_value].iterator $it = $input_cpp_var.begin()
                |cdef $cy_tt $item
                |while $it != $input_cpp_var.end():
                |   $item = $cy_tt.__new__($cy_tt)
                |   $item.inst = shared_ptr[$cy_tt_value](new $cy_tt_value((deref($it)).second))
                |   $output_py_var[$key_conv] = $item
                |   inc($it)
                """, locals())
            return code
        else:
            value_conv = "<%s>(deref(%s).second)" % (cy_tt_value, it)
            code = Code().add("""
                |$output_py_var = dict()
                |cdef _Map[$cy_tt_key, $cy_tt_value].iterator $it = $input_cpp_var.begin()
                |while $it != $input_cpp_var.end():
                |   $output_py_var[$key_conv] = $value_conv
                |   inc($it)
                """, locals())
            return code


import time

class CVTermMapConverter(TypeConverterBase):

    def get_base_types(self):
        return "Map",

    def matches(self, cpp_type):
        print str(cpp_type), "Map[String,libcpp_vector[CVTerm]]"
        return str(cpp_type) == "Map[String,libcpp_vector[CVTerm]]" \
           or  str(cpp_type) == "Map[String,libcpp_vector[CVTerm]] &"

    def matching_python_type(self, cpp_type):
        return "dict"

    def type_check_expression(self, cpp_type, arg_var):
        return Code().add("""
          |isinstance($arg_var, dict)
          + and all(isinstance(k, str) for k in $arg_var.keys())
          + and all(isinstance(v, list) for v in $arg_var.values())
          + and all(isinstance(vi, CVTerm) for v in $arg_var.values() for vi in
          v)
          """, locals()).render()

    def input_conversion(self, cpp_type, argument_var, arg_num):

        map_name = "_map_%d" % arg_num
        v_vec = "_v_vec_%d" % arg_num
        v_ptr = "_v_ptr_%d" % arg_num
        v_i = "_v_i_%d" % arg_num
        k_string = "_k_str_%d" % arg_num

        code = Code().add("""
                |cdef Map[_String, libcpp_vector[_CVTerm]] $map_name
                |cdef libcpp_vector[_CVTerm] $v_vec
                |cdef _String $k_string
                |cdef CVTerm $v_i
                |for k, v in $argument_var.items():
                |    $v_vec.clear()
                |    for $v_i in v:
                |        $v_vec.push_back(deref($v_i.inst.get()))
                |    $map_name[_String(<char *>k)] = $v_vec
                """, locals())

        if cpp_type.is_ref:
            replace = "_replace_%d" % arg_num
            outer_it = "outer_it_%d" % arg_num
            inner_it = "inner_it_%d" % arg_num
            item     = "item_%d" % arg_num
            inner_key = "inner_key_%d" % arg_num
            inner_values = "inner_values_%d" % arg_num
            cleanup_code = Code().add("""
                |cdef $replace = dict()
                |cdef Map[_String, libcpp_vector[_CVTerm]].iterator $outer_it
                + = $map_name.begin()

                |cdef libcpp_vector[_CVTerm].iterator $inner_it
                |cdef CVTerm $item
                |cdef str $inner_key
                |cdef list $inner_values

                |while $outer_it != $map_name.end():
                |   $inner_key = deref($outer_it).first.c_str()
                |   $inner_values = []
                |   $inner_it = deref($outer_it).second.begin()
                |   while $inner_it != deref($outer_it).second.end():
                |       $item = CVTerm.__new__(CVTerm)
                |       $item.inst = shared_ptr[_CVTerm](new _CVTerm(deref($inner_it)))
                |       $inner_values.append($item)
                |       inc($inner_it)
                |   $replace[$inner_key] = $inner_values
                |   inc($outer_it)

                |$argument_var.clear()
                |$argument_var.update($replace)
                """, locals())
        else:
            cleanup_code = ""
        return code, map_name, cleanup_code


    def call_method(self, res_type, cy_call_str):
        return "_r = %s" % (cy_call_str)

    def output_conversion(self, cpp_type, input_cpp_var, output_py_var):

        rnd = hex(id(self))+str(time.time()).split(".")[0]
        outer_it = "outer_it_%s" % rnd
        inner_it = "inner_it_%s" % rnd
        item     = "item_%s" % rnd
        inner_key = "inner_key_%s" % rnd
        inner_values = "inner_values_%s" % rnd

        code = Code().add("""
            |$output_py_var = dict()
            |cdef Map[_String, libcpp_vector[_CVTerm]].iterator $outer_it
            + = $input_cpp_var.begin()
            |cdef libcpp_vector[_CVTerm].iterator $inner_it
            |cdef CVTerm $item
            |cdef str $inner_key
            |cdef list $inner_values
            |while $outer_it != $input_cpp_var.end():
            |   $inner_key = deref($outer_it).first.c_str()
            |   $inner_values = []
            |   $inner_it = deref($outer_it).second.begin()
            |   while $inner_it != deref($outer_it).second.end():
            |       $item = CVTerm.__new__(CVTerm)
            |       $item.inst = shared_ptr[_CVTerm](new _CVTerm(deref($inner_it)))
            |       $inner_values.append($item)
            |       inc($inner_it)
            |   $output_py_var[$inner_key] = $inner_values
            |   inc($outer_it)
            """, locals())

        return code


