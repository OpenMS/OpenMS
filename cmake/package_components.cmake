cpack_add_install_type(recommended DISPLAY_NAME "Recommended")
cpack_add_install_type(full DISPLAY_NAME "Full")
cpack_add_install_type(minimal DISPLAY_NAME "Minimal")

cpack_add_component(share
                DISPLAY_NAME "OpenMS shared files"
                DESCRIPTION "OpenMS shared files"
                INSTALL_TYPES recommended full minimal
                )
cpack_add_component(library
                DISPLAY_NAME "Libraries"
                DESCRIPTION "Libraries"
                INSTALL_TYPES recommended full minimal
                )
cpack_add_component(applications
                DISPLAY_NAME "OpenMS binaries"
                DESCRIPTION "OpenMS binaries including TOPP tools/utils, TOPPView and TOPPAS."
                INSTALL_TYPES recommended full minimal
                )
cpack_add_component(doc
                DISPLAY_NAME "Documentation"
                DESCRIPTION "Class and tool documentation. With tutorials."
                INSTALL_TYPES recommended full
                )
cpack_add_component_group(thirdparty
                     DISPLAY_NAME "Thirdparty binaries"
                     DESCRIPTION "Binaries and files for thirdparty tools and engines."
                     EXPANDED
                     )
foreach(component IN LISTS ${THIRDPARTY_COMPONENT_GROUP})
    cpack_add_component(${component}
                    DISPLAY_NAME ${component}
                    DESCRIPTION "Thirdparty engine ${component}"
                    GROUP thirdparty
                    INSTALL_TYPES recommended full
                    )
endforeach()
