#!/bin/bash

# 2012/01/06 - Stephan Aiche
# based on http://stackoverflow.com/questions/96882/how-do-i-create-a-nice-looking-dmg-for-mac-os-x-using-command-line-tools
# bash script to further customize the OpenMS dmg

# if the apple script part fails with "timed out" error try
# defaults write com.apple.Finder AppleShowAllFiles TRUE && killall Finder
# <wait 1-2min>
# defaults write com.apple.Finder AppleShowAllFiles FALSE && killall Finder

DISK_NAME=OpenMS-2.3.0-Darwin
DMG_NAME=${DISK_NAME}.dmg
OPENMS_NAME=OpenMS-2.3.0
backgroundPictureName=.background.png
LICENSE=${PWD}/_CPack_Packages/Darwin/DragNDrop/sla.r

if [ ! -e ${DMG_NAME} ]
then
  echo "Please execute make package first"
  exit 1
fi

# make dmg writable
hdiutil convert ${DMG_NAME} -format UDRW -o temp.dmg
#hdiutil attach temp.dmg

device=$(hdiutil attach -readwrite -noverify -noautoopen "temp.dmg" | \
         egrep '^/dev/' | sed 1q | awk '{print $1}')

# remove original dmg
rm -f ${DMG_NAME}

# wait till package is open
#sleep 10

echo 'tell application "Finder"
	tell disk "'${DISK_NAME}'"
    with timeout of 300 seconds
  		open

			set theXOrigin to 400
			set theYOrigin to 200
			set theBottomRightX to 1030
			set theBottomRightY to 785

			tell container window
				set current view to icon view
				set toolbar visible to false
				set statusbar visible to false
				set the bounds to {theXOrigin, theYOrigin, theBottomRightX, theBottomRightY}
				set statusbar visible to false
			end tell

			update without registering applications
			delay 1

  		set theViewOptions to the icon view options of container window
  		set arrangement of theViewOptions to not arranged

      -- if we have a fixed resolution we can also set this
    	set icon size of theViewOptions to 72
  		set bgimg to "'${OPENMS_NAME}':share:OpenMS:background.png" as text
  		set background picture of theViewOptions to file bgimg

  		set the position of item "'${OPENMS_NAME}'" of container window to {470, 140}
  		set the position of item "Applications" of container window to {160, 140}

  		--give the finder some time to write the .DS_Store file
  		delay 3

  		-- work around for Snow Leopard bug
  		close
  		open

  		update without registering applications
  		-- delay 5

      -- try eject
      eject
    end timeout
	end tell
end tell
' | osascript

# be sure that everything is done
#sleep 10

chmod -Rf go-w /Volumes/"${DISK_NAME}"
sync
sync
hdiutil detach ${device}
hdiutil convert "temp.dmg" -format UDZO -imagekey zlib-level=9 -o "${DMG_NAME}"
#rm -f /pack.temp.dmg

# hdiutil convert temp.dmg -format UDZO -imagekey zlib-level=9 -o ${DMG_NAME}
rm -f temp.dmg

if [ ! -z "${LICENSE}" -a "${LICENSE}" != "-null-" ]; then
  echo "adding EULA resources"
  hdiutil unflatten "${DMG_NAME}"
  ResMerger "${LICENSE}" -a -o "${DMG_NAME}"
  hdiutil flatten "${DMG_NAME}"
fi
