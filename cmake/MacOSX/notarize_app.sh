#!/bin/bash
set -e

# Modern notarization script using notarytool (Xcode 13+)
# Usage: notarize_app.sh <bundle_path> <bundle_id> <apple_id> <password_env_var> <log_folder>

BUNDLE_PKG="$1"
BUNDLE_ID="$2"
APPLE_ID="$3"
PASSWORD_ENV_VAR="$4"
LOG_FOLDER="$5"

if [ -z "$BUNDLE_PKG" ] || [ -z "$BUNDLE_ID" ] || [ -z "$APPLE_ID" ] || [ -z "$PASSWORD_ENV_VAR" ]; then
  echo "Error: Missing required arguments"
  echo "Usage: $0 <bundle_path> <bundle_id> <apple_id> <password_env_var> <log_folder>"
  exit 1
fi

# Get the password from the environment variable
APP_PASSWORD="${!PASSWORD_ENV_VAR}"
if [ -z "$APP_PASSWORD" ]; then
  echo "Error: Password environment variable $PASSWORD_ENV_VAR is not set"
  exit 1
fi

NOTARIZE_LOG="${LOG_FOLDER}/notarize.log"
mkdir -p "$LOG_FOLDER"

REMOVE_PKG=false
BUNDLE_FILE="$BUNDLE_PKG"

# Handle different file types
if [[ $BUNDLE_PKG == *.dmg ]] || [[ $BUNDLE_PKG == *.pkg ]]; then
  # dmg and pkg can be notarized directly
  BUNDLE_FILE="$BUNDLE_PKG"
elif [[ $BUNDLE_PKG == *.zip ]]; then
  # zip files: assume the bundle is in the same folder for stapling
  BUNDLE_FILE="${BUNDLE_PKG%.*}"
elif [[ $BUNDLE_PKG == *.app ]]; then
  # app bundles: need to be zipped for upload
  BUNDLE_FILE="$BUNDLE_PKG"
  BUNDLE_PKG="${BUNDLE_PKG}.zip"
  echo "Creating zip archive for app bundle..."
  ditto -c -k --rsrc --keepParent "$BUNDLE_FILE" "$BUNDLE_PKG"
  REMOVE_PKG=true
else
  echo "Error: Unsupported file type for notarization: $BUNDLE_PKG"
  exit 1
fi

echo "Submitting $BUNDLE_PKG for notarization..."

# Submit for notarization using notarytool (modern approach)
# The --wait flag makes notarytool wait for the notarization to complete
if xcrun notarytool submit "$BUNDLE_PKG" \
  --apple-id "$APPLE_ID" \
  --password "$APP_PASSWORD" \
  --team-id "C64UCGJ5PL" \
  --wait \
  --timeout 30m \
  2>&1 | tee "$NOTARIZE_LOG"; then
  
  echo "Notarization successful!"
  
  # Staple the notarization ticket to the bundle
  echo "Stapling notarization ticket to $BUNDLE_FILE..."
  if xcrun stapler staple "$BUNDLE_FILE" 2>&1 | tee -a "$NOTARIZE_LOG"; then
    echo "Stapling successful!"
  else
    echo "Warning: Stapling failed, but notarization was successful"
    # Note: Stapling can fail for some file types (e.g., zip files)
    # but the notarization ticket is still available online
  fi
  
  # Clean up temporary zip if created
  if [ "$REMOVE_PKG" = true ]; then
    rm "$BUNDLE_PKG"
  fi
  
  exit 0
else
  echo "Error: Notarization failed"
  cat "$NOTARIZE_LOG"
  
  # Clean up temporary zip if created
  if [ "$REMOVE_PKG" = true ]; then
    rm "$BUNDLE_PKG"
  fi
  
  exit 1
fi
