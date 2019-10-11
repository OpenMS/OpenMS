#!/bin/bash
#set -x

ASC_USERNAME="$3"
ASC_PASSWORD_ENVVAR="$4"

BUNDLE_ID="$2"
BUNDLE_PKG="$1"

LOG_FOLDER="$5"

NOTARIZE_INFO_LOG=$LOG_FOLDER/notarize_info.log
NOTARIZE_APP_LOG=$LOG_FOLDER/notarize_app.log

touch $NOTARIZE_INFO_LOG
touch $NOTARIZE_APP_LOG

# submit app for notarization
if xcrun altool --notarize-app --primary-bundle-id "$BUNDLE_ID" --username "$ASC_USERNAME" --password @env:$ASC_PASSWORD_ENVVAR -f "$BUNDLE_PKG" > "$NOTARIZE_APP_LOG" 2>&1; then
	cat "$NOTARIZE_APP_LOG"
	RequestUUID=$(awk -F ' = ' '/RequestUUID/ {print $2}' "$NOTARIZE_APP_LOG")

	# check status periodically
	while sleep 60 && date; do
		# check notarization status
		if xcrun altool --notarization-info "$RequestUUID" --username "$ASC_USERNAME" --password @env:$ASC_PASSWORD_ENVVAR > "$NOTARIZE_INFO_LOG" 2>&1; then
			cat "$NOTARIZE_INFO_LOG"

			# once notarization is complete, run stapler and exit
			if ! grep -q "Status: in progress" "$NOTARIZE_INFO_LOG"; then
				xcrun stapler staple "$BUNDLE_PKG"
				exit
			fi
		else
			cat "$NOTARIZE_INFO_LOG" 1>&2
			exit 1
		fi
	done
else
	cat "$NOTARIZE_APP_LOG" 1>&2
	exit 1
fi
