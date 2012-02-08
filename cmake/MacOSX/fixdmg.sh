#!/bin/bash

# 2012/01/06 - Stephan Aiche
# based on http://stackoverflow.com/questions/96882/how-do-i-create-a-nice-looking-dmg-for-mac-os-x-using-command-line-tools
# bash script to further customize the OpenMS dmg

DISK_NAME=OpenMS-1.9.0-Darwin
DMG_NAME=${DISK_NAME}.dmg
OPENMS_NAME=OpenMS-1.9.0
backgroundPictureName=.background.png

if [ ! -e ${DMG_NAME} ]
then
  echo "Please execute make package first"
  exit 1
fi

# make dmg writable
hdiutil convert ${DMG_NAME} -format UDRW -o temp.dmg
hdiutil attach temp.dmg

# remove original dmg
rm -f ${DMG_NAME}

# wait till package is open
sleep 10

echo '
  tell application "Finder"
	tell disk "'${DISK_NAME}'"
		open
		set current view of container window to icon view
		set toolbar visible of container window to false
		set statusbar visible of container window to false
		set the bounds of container window to {400, 200, 1030, 785}
		set theViewOptions to the icon view options of container window
		set arrangement of theViewOptions to not arranged
    -- if we have a fixed resolution we can also set this
  	set icon size of theViewOptions to 72
		set bgimg to "'${OPENMS_NAME}':share:OpenMS:background.png" as text
		set background picture of theViewOptions to file bgimg
		
		set the position of item "'${OPENMS_NAME}'" of container window to {470, 140}
		set the position of item "Applications" of container window to {160, 140}
		
		-- work around for Snow Leopard bug
		close
		open
		
		update without registering applications
		delay 5
		eject
	end tell
end tell
' | osascript

# be sure that everything is done
sleep 10

hdiutil convert temp.dmg -format UDZO -imagekey zlib-level=9 -o ${DMG_NAME}
rm -f temp.dmg
