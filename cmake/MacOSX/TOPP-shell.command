#!/bin/bash
clear
echo ""
echo ""
echo "   ------------------------------------------------"
echo "       Welcome to the OpenMS/TOPP command shell    "
echo "   ------------------------------------------------"
echo ""
echo ""

export OPENMS_TOPP_PATH=`dirname "$0"`
# declare as TOPP_SHELL
export IS_TOPP_SHELL="TOPP-SHELL"
bash --rcfile ${OPENMS_TOPP_PATH}/.TOPP_bash_profile

