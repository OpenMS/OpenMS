



    def encode64(self, list in_ , int to_byte_order ,  out ,  zlib_compression ):
        """Cython signature: void encode64(libcpp_vector[double] & in_, ByteOrder to_byte_order, String & out, bool zlib_compression)"""
        assert isinstance(in_, list) and all(isinstance(elemt_rec, float) for elemt_rec in in_), 'arg in_ wrong type'
        assert to_byte_order in [0, 1], 'arg to_byte_order wrong type'
        assert isinstance(out, String), 'arg out wrong type'
        assert isinstance(zlib_compression, (int, long)), 'arg zlib_compression wrong type'
        cdef libcpp_vector[double] v0 = in_
    
    
    
        self.inst.get().encode(v0, (<_ByteOrder>to_byte_order), deref((<String>out).inst.get()), (<bool>zlib_compression))
    
    def decode64(self,  in_ , int from_byte_order , list out ,  zlib_compression ):
        """Cython signature: void decode64(const String & in_, ByteOrder from_byte_order, libcpp_vector[double] & out, bool zlib_compression)"""
        assert (isinstance(in_, str) or isinstance(in_, unicode) or isinstance(in_, bytes) or isinstance(in_, String)), 'arg in_ wrong type'
        assert from_byte_order in [0, 1], 'arg from_byte_order wrong type'
        assert isinstance(out, list) and all(isinstance(elemt_rec, float) for elemt_rec in out), 'arg out wrong type'
        assert isinstance(zlib_compression, (int, long)), 'arg zlib_compression wrong type'
    
    
        cdef libcpp_vector[double] v2 = out
    
        self.inst.get().decode(deref((convString(in_)).get()), (<_ByteOrder>from_byte_order), v2, (<bool>zlib_compression))
        out[:] = v2

    def encode32(self, list in_ , int to_byte_order ,  out ,  zlib_compression ):
        """Cython signature: void encode32(libcpp_vector[float] & in_, ByteOrder to_byte_order, String & out, bool zlib_compression)"""
        assert isinstance(in_, list) and all(isinstance(elemt_rec, float) for elemt_rec in in_), 'arg in_ wrong type'
        assert to_byte_order in [0, 1], 'arg to_byte_order wrong type'
        assert isinstance(out, String), 'arg out wrong type'
        assert isinstance(zlib_compression, (int, long)), 'arg zlib_compression wrong type'
        cdef libcpp_vector[float] v0 = in_
    
    
    
        self.inst.get().encode(v0, (<_ByteOrder>to_byte_order), deref((<String>out).inst.get()), (<bool>zlib_compression))
    
    def decode32(self,  in_ , int from_byte_order , list out ,  zlib_compression ):
        """Cython signature: void decode32(const String & in_, ByteOrder from_byte_order, libcpp_vector[float] & out, bool zlib_compression)"""
        assert (isinstance(in_, str) or isinstance(in_, unicode) or isinstance(in_, bytes) or isinstance(in_, String)), 'arg in_ wrong type'
        assert from_byte_order in [0, 1], 'arg from_byte_order wrong type'
        assert isinstance(out, list) and all(isinstance(elemt_rec, float) for elemt_rec in out), 'arg out wrong type'
        assert isinstance(zlib_compression, (int, long)), 'arg zlib_compression wrong type'
    
    
        cdef libcpp_vector[float] v2 = out
    
        self.inst.get().decode(deref((convString(in_)).get()), (<_ByteOrder>from_byte_order), v2, (<bool>zlib_compression))
        out[:] = v2
