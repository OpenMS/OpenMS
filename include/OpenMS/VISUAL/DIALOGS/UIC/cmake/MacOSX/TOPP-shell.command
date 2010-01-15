#!/bin/bash
clear
echo ""
echo ""
echo "   ------------------------------------------------"
echo "       Welcome to the OpenMS/TOPP command shell    "
echo "   ------------------------------------------------"
echo ""
echo ""

export OPENMS_TOPP_PATH=`dirname $0`
bash --rcfile ${OPENMS_TOPP_PATH}/.TOPP_bash_profile

