#!/bin/bash
# macOS notarization script using notarytool (requires Xcode 13+ / macOS 11.3+)
# altool was deprecated and unsupported after Fall 2023.
#
# Usage: notarize.sh <bundle_pkg> <bundle_id> <apple_id> <password_env_var> [log_folder]
#
# Arguments:
#   bundle_pkg       - The package to notarize (.dmg, .pkg, .zip, or .app)
#   bundle_id        - The bundle identifier (e.g., de.openms)
#   apple_id         - Apple ID email for notarization
#   password_env_var - Environment variable name containing app-specific password
#   log_folder       - Optional: folder for log files (defaults to current directory)

# Exit on error and fail on any error in a pipeline
set -e
set -o pipefail

BUNDLE_PKG="$1"
BUNDLE_ID="$2"
ASC_USERNAME="$3"
ASC_PASSWORD_ENVVAR="$4"
ASC_TEAMID="$5"
LOG_FOLDER="${6:-.}"

NOTARIZE_LOG="$LOG_FOLDER/notarize.log"

mkdir -p "$LOG_FOLDER"
touch "$NOTARIZE_LOG"

REMOVE_PKG=false
IS_ZIP=false

echo "=== macOS Notarization Script ==="
echo "Bundle: $BUNDLE_PKG"
echo "Bundle ID: $BUNDLE_ID"
echo "Apple ID: $ASC_USERNAME"
echo "Log folder: $LOG_FOLDER"

# Validate inputs
if [[ -z "$BUNDLE_PKG" ]] || [[ -z "$BUNDLE_ID" ]] || [[ -z "$ASC_USERNAME" ]] || [[ -z "$ASC_PASSWORD_ENVVAR" ]]; then
    echo "Error: Missing required arguments"
    echo "Usage: $0 <bundle_pkg> <bundle_id> <apple_id> <password_env_var> [log_folder]"
    exit 1
fi

# Check that the password environment variable is set
if [[ -z "${!ASC_PASSWORD_ENVVAR}" ]]; then
    echo "Error: Environment variable '$ASC_PASSWORD_ENVVAR' is not set or empty"
    exit 1
fi

# Determine file type and prepare for notarization
# Only notarize the outermost container (zip, pkg, dmg)
if [[ $BUNDLE_PKG == *.dmg ]]; then
    BUNDLE_FILE=$BUNDLE_PKG
    echo "Notarizing DMG: $BUNDLE_PKG"
elif [[ $BUNDLE_PKG == *.pkg ]]; then
    BUNDLE_FILE=$BUNDLE_PKG
    echo "Notarizing PKG: $BUNDLE_PKG"
elif [[ $BUNDLE_PKG == *.zip ]]; then
    # For zip files, we need to unzip to staple, then re-zip
    BUNDLE_FILE="$BUNDLE_PKG"
    IS_ZIP=true
    echo "Notarizing ZIP: $BUNDLE_PKG (will staple contents)"
elif [[ $BUNDLE_PKG == *.app ]]; then
    # Apps need to be zipped for upload, then unzipped for stapling
    BUNDLE_FILE=$BUNDLE_PKG
    BUNDLE_PKG="${BUNDLE_PKG}.zip"
    echo "Zipping app bundle for notarization..."
    ditto -c -k --rsrc --keepParent "$BUNDLE_FILE" "$BUNDLE_PKG"
    REMOVE_PKG=true
    echo "Notarizing APP (via ZIP): $BUNDLE_PKG"
else
    echo "Error: Unsupported filetype for notarization: $BUNDLE_PKG"
    echo "Supported types: .dmg, .pkg, .zip, .app"
    exit 1
fi

# Verify the file exists
if [[ ! -f "$BUNDLE_PKG" ]]; then
    echo "Error: File not found: $BUNDLE_PKG"
    exit 1
fi

echo ""
echo "=== Submitting for notarization ==="

# Submit for notarization using notarytool
# --wait makes the command block until notarization is complete
# Ensure pipefail is set for this block in case of subshells
set -o pipefail
if xcrun notarytool submit "$BUNDLE_PKG" \
    --apple-id "$ASC_USERNAME" \
    --password "${!ASC_PASSWORD_ENVVAR}" \
    --team-id "$ASC_TEAMID" \
    --wait \
    2>&1 | tee "$NOTARIZE_LOG"; then

    echo ""
    echo "=== Notarization submission completed ==="

    # Check if notarization was successful
    if grep -q "status: Accepted" "$NOTARIZE_LOG"; then
        echo "Notarization successful!"

        # Staple the notarization ticket to the bundle
        echo ""
        echo "=== Stapling notarization ticket ==="

        # Note: You cannot staple a .zip file directly
        # If the original was a zip, we need to handle it differently
        if [[ "$IS_ZIP" = true ]]; then
            echo "Warning: Cannot staple a .zip file. The notarization is stored with Apple."
            echo "Users will need to be online for Gatekeeper to verify the notarization."
        else
            if xcrun stapler staple "$BUNDLE_FILE"; then
                echo "Stapling successful!"

                # Verify the stapling worked
                echo ""
                echo "=== Verifying notarization ==="
                xcrun stapler validate "$BUNDLE_FILE" || echo "Warning: Stapler validation returned non-zero"
            else
                echo "Warning: Stapling failed, but notarization was successful."
                echo "Users will need to be online for Gatekeeper to verify."
            fi
        fi

        # Clean up temporary zip if we created one
        if [ "$REMOVE_PKG" = true ]; then
            rm -f "$BUNDLE_PKG"
        fi

        echo ""
        echo "=== Notarization complete ==="
        exit 0
    else
        echo "Error: Notarization failed!"
        echo "Check the log for details: $NOTARIZE_LOG"

        # Try to get the submission ID and fetch detailed logs
        SUBMISSION_ID=$(grep -o 'id: [a-f0-9-]*' "$NOTARIZE_LOG" | head -1 | cut -d' ' -f2)
        if [[ -n "$SUBMISSION_ID" ]]; then
            echo ""
            echo "=== Fetching detailed notarization log ==="
            xcrun notarytool log "$SUBMISSION_ID" \
                --apple-id "$ASC_USERNAME" \
                --password "${!ASC_PASSWORD_ENVVAR}" \
                "$LOG_FOLDER/notarization_details.json" 2>&1 || true

            if [[ -f "$LOG_FOLDER/notarization_details.json" ]]; then
                echo "Detailed log saved to: $LOG_FOLDER/notarization_details.json"
                cat "$LOG_FOLDER/notarization_details.json"
            fi
        fi

        # Clean up temporary zip if we created one
        if [ "$REMOVE_PKG" = true ]; then
            rm -f "$BUNDLE_PKG"
        fi

        exit 1
    fi
else
    echo "Error: notarytool submission failed!"
    cat "$NOTARIZE_LOG"

    # Clean up temporary zip if we created one
    if [ "$REMOVE_PKG" = true ]; then
        rm -f "$BUNDLE_PKG"
    fi

    exit 1
fi
