#DEPRECATED. Functionality has been merged into update_version_numbers.yml
TOOL_DIR_PATH="./src/tests/topp/"
find $TOOL_DIR_PATH -type f -iname '*.ini' -exec grep -q '<ITEM name="version".*Version of the tool' {} \; -exec sed -i '' -e 's/name="version" value=".*" type="string"/name="version" value="3.2.0" type="string"/g' {} \;
