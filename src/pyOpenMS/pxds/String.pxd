from libc.string cimport const_char
from libcpp cimport bool
from libcpp.string cimport string as libcpp_string

cdef extern from "<OpenMS/DATASTRUCTURES/String.h>" namespace "OpenMS":

    # Note that all input will be encoded as UTF8 and stored as bytes inside
    # OpenMS::String (which is a wrapper around std::string). The data will be
    # decoded on its way back when calling "toString()" again as UTF8.
    # If you use another encoding, you will have a bad time!
    # Also, note that __str__ follows Python 2.x convention and returns bytes,
    # not unicode. Please use toString if you want to have encoded data.
    cdef cppclass String:
        # wrap-hash:
        #   c_str()

        String() nogil except +
        String(String) nogil except +  # wrap-ignore
        String(char *) nogil except + # wrap-ignore
        String(char *, size_t l) nogil except + # wrap-ignore
        String(str) nogil except + # wrap-ignore
        const_char * c_str() nogil except + # wrap-ignore

        # Creates a Python 2/3 unicode string (use this instead of c_str() if you
        # plan to use any non-ASCII code).
        toString(self) nogil except + # wrap-ignore

        bool operator==(String) nogil except +
        bool operator!=(String) nogil except +

        # Rather perform string operations in Python (you will have a bad time
        # with unicode strings otherwise).
        size_t length() nogil except + # wrap-ignore
        # libcpp_string operator[](int) nogil except + # wrap-upper-limit:length()

cdef extern from "<OpenMS/DATASTRUCTURES/String.h>" namespace "OpenMS::String":
    
    cdef enum QuotingMethod "OpenMS::String::QuotingMethod":
        NONE
        ESCAPE
        DOUBLE


#
#   def _init_2(self, str in_0 ):
#       assert isinstance(in_0, str), 'arg in_0 wrong type'
#   
#       print(b"have init 2 here")
#       # TODO catch encoding errors ...
#       py_byte_string = in_0.encode('UTF-8')
#       cdef char* c_string = py_byte_string
#       self.inst = shared_ptr[_String](new _String((<char *>c_string)))
#
#       # self.inst = shared_ptr[_String](new _String((<char *>in_0)))
#   
#   def __init__(self, *args):
#       if not args:
#            self._init_0(*args)
#       elif (len(args)==1) and (isinstance(args[0], bytes)):
#            self._init_1(*args)
#       elif (len(args)==1) and (isinstance(args[0], str)):
#            self._init_2(*args)
#       else:
#              print(b"aaaa")
#              raise Exception('can not handle type of %s' % (args,)) 
#
