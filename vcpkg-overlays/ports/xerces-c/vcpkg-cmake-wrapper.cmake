# FindXercesC omits transitive dependencies such as curl from static libraries.
list(REMOVE_ITEM ARGS "NO_MODULE" "CONFIG" "MODULE")
_find_package(${ARGS} CONFIG)

if (APPLE)
   if (TARGET XercesC::XercesC)
      set_property(TARGET XercesC::XercesC APPEND PROPERTY INTERFACE_LINK_LIBRARIES  "-framework CoreServices" "-framework CoreFoundation")
      list(APPEND XercesC_LIBRARIES "-framework CoreServices" "-framework CoreFoundation")
   endif()
endif()
