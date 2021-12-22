#!/bin/bash
#set -x

ASC_USERNAME="$3"
## currently only supports providing PW through an env. variable
ASC_PASSWORD_ENVVAR="$4"

BUNDLE_ID="$2"
BUNDLE_PKG="$1"

LOG_FOLDER="$5"

NOTARIZE_INFO_LOG=$LOG_FOLDER/notarize_info.log
NOTARIZE_APP_LOG=$LOG_FOLDER/notarize_app.log

touch $NOTARIZE_INFO_LOG
touch $NOTARIZE_APP_LOG

REMOVE_PKG=false

## Checks if it is a zip, then the unzipped bundle should be in the same folder for stapling
## if it is e.g. a dmg already, nothing needs to be done
if [[ $BUNDLE_PKG == *dmg ]] ; then
  BUNDLE_FILE=$BUNDLE_PKG
elif [[ $BUNDLE_PKG == *zip ]] ; then
  ## assume file is in same folder
  BUNDLE_FILE=${BUNDLE_PKG%.*}
elif [[ $BUNDLE_PKG == *app ]] ; then
  ## zip it for upload and remove later
  BUNDLE_FILE=$BUNDLE_PKG
  BUNDLE_PKG=${BUNDLE_PKG}.zip
  ditto -c -k --rsrc --keepParent $BUNDLE_FILE $BUNDLE_PKG
  REMOVE_PKG=true
else
  echo "Unsupported filetype for notarization. $BUNDLE_PKG"
  exit 1
fi



# submit app for notarization
if xcrun altool --notarize-app --primary-bundle-id "$BUNDLE_ID" --username "$ASC_USERNAME" --password @env:$ASC_PASSWORD_ENVVAR -f "$BUNDLE_PKG" > "$NOTARIZE_APP_LOG" 2>&1; then
	cat "$NOTARIZE_APP_LOG"
	RequestUUID=$(awk -F ' = ' '/RequestUUID/ {print $2}' "$NOTARIZE_APP_LOG")

	# check status periodically
	while sleep 60 && date; do
		# check notarization status
		if xcrun altool --notarization-info "$RequestUUID" --username "$ASC_USERNAME" --password @env:$ASC_PASSWORD_ENVVAR > "$NOTARIZE_INFO_LOG" 2>&1; then
			cat "$NOTARIZE_INFO_LOG"

			# once notarization is complete, run stapler on the unzipped bundle and exit
			if ! grep -q "Status: in progress" "$NOTARIZE_INFO_LOG"; then
				xcrun stapler staple "$BUNDLE_FILE"
				if [ "$REMOVE_PKG" = true ] ; then
				  rm $BUNDLE_PKG
				fi
				exit
			fi
		else
			cat "$NOTARIZE_INFO_LOG" 1>&2
			if [ "$REMOVE_PKG" = true ] ; then
			  rm $BUNDLE_PKG
			fi
			exit 1
		fi
	done
else
	cat "$NOTARIZE_APP_LOG" 1>&2
	if [ "$REMOVE_PKG" = true ] ; then
	  rm $BUNDLE_PKG
	fi	
	exit 1
fi
