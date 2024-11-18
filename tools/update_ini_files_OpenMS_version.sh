TOOL_DIR_PATH="./src/tests/topp/"
OPENMS_VERSION="3.2.0"
find $TOOL_DIR_PATH -type f -iname '*.ini' \
  -exec grep -q '<ITEM name="version".*Version of the tool' {} \; \
  -exec sed -i -e "s/name=\"version\" value=\".*\" type=\"string\"/name=\"version\" value=\"$OPENMS_VERSION\" type=\"string\"/g" {} \;
