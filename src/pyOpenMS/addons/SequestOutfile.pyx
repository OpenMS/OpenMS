



    def getSequences(self, String database_filename , dict ac_position_map, list sequences, list found, dict not_found):
        # void getSequences(String &database_filename,
        # libcpp_map[ String, Size ] &ac_position_map,
        # libcpp_vector[ String ] &sequences,
        # libcpp_vector[ libcpp_pair[ String, Size ] ] &found,
        # libcpp_map[ String, Size ] &not_found) nogil except +
        assert isinstance(database_filename, String), 'arg database_filename wrong type'
        assert isinstance(sequences, list) and all(isinstance(i, bytes) for i in sequences), 'arg sequences wrong type'
    
        cdef libcpp_vector[_String] * v1 = new libcpp_vector[_String]()
        cdef bytes item1
        for item1 in sequences:
           v1.push_back(_String(<char *>item1))

        cdef libcpp_map[_String, size_t ] c_dict_ac_position
        for k,v in ac_position_map.iteritems():
            c_dict_ac_position[ _String(<char *>k) ] = v

        cdef libcpp_map[_String, size_t ] c_dict_not_found
        for k,v in not_found.iteritems():
            c_dict_not_found[ _String(<char *>k) ] = v

        cdef libcpp_vector[ libcpp_pair[_String, size_t] ] * v4 = new libcpp_vector[ libcpp_pair[_String, size_t] ]()
        cdef list item4
        cdef libcpp_pair[_String, size_t] * aPair
        for item4 in found:

           aPair = new libcpp_pair[ _String, size_t] ( _String(<char *>item4[0]), item4[1])
           v4.push_back( deref(aPair) )
           del aPair
 
        #
        # Call C++
        #
        self.inst.get().getSequences(deref(database_filename.inst.get()), c_dict_ac_position, deref(v1), deref(v4), c_dict_not_found)
        #
        # Get arguments back from C++
        #

        # libcpp_map[ String, Size ] &ac_position_map,
        replace_2 = dict()
        cdef libcpp_map[_String, size_t].iterator it_dict2 = c_dict_ac_position.begin()
        while it_dict2 != c_dict_ac_position.end():
            replace_2[ <libcpp_string>deref(it_dict2).first ] = deref(it_dict2).second
            inc(it_dict2)
        ac_position_map.update(replace_2)

        # libcpp_vector[ String ] &sequences,
        cdef replace = []
        cdef libcpp_vector[_String].iterator it = v1.begin()
        while it != v1.end():
           replace.append(<char*>deref(it).c_str())
           inc(it)
        sequences[:] = replace

        # libcpp_vector[ libcpp_pair[ String, Size ] ] &found,
        replace_4 = list()
        cdef libcpp_vector[ libcpp_pair[_String, size_t] ].iterator it_list4 = deref(v4).begin()
        cdef libcpp_pair[_String, size_t] anotherPair
        while it_list4 != deref(v4).end():
            # get the pair back
            anotherPair = deref(it_list4)
            replace_4.append( [<libcpp_string>anotherPair.first, anotherPair.second] )
            inc(it_list4)
        found[:] = replace_4

        # libcpp_map[ String, Size ] &not_found
        replace_5 = dict()
        cdef libcpp_map[_String, size_t].iterator it_dict5 = c_dict_not_found.begin()
        while it_dict5 != c_dict_not_found.end():
            replace_2[ <libcpp_string>deref(it_dict5).first ] = deref(it_dict5).second
            inc(it_dict5)
        not_found.update(replace_2)


        del v1
        del v4

