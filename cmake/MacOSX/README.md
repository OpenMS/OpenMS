### Installing TOPP into your local shell ###

To use TOPP tools as regular applications in your own shell just add the following lines to your ~/.profile file. 
(Adapt the first line to point to the folder where you put OpenMS, e.g.: /Applications/OpenMS-3.0.0)

```shell
export OPENMS_TOPP_PATH=<OpenMS-PATH>
source $OPENMS_TOPP_PATH/.TOPP_bash_profile
```

---

## macOS Code Signing and Notarization

OpenMS packages for macOS are code signed and notarized to ensure users can run them without security warnings from Gatekeeper.

### Prerequisites

1. **Apple Developer Account** ($99/year) - https://developer.apple.com
2. **Xcode** installed with developer tools (`xcode-select --install`)
3. **Certificates** exported from Xcode:
   - **Developer ID Application** - for signing binaries
   - **Developer ID Installer** - for signing .pkg installers

### Required GitHub Secrets

For CI/CD, the following secrets must be configured in your GitHub repository.

**Note:** Only certificates and passwords need to be secrets. Identity strings and email addresses are public information.

| Secret | Description |
|--------|-------------|
| `APPLE_DEVELOPER_CERT_64` | Base64-encoded Developer ID Application .p12 certificate |
| `APPLE_DEVELOPMENT_PW` | Password for the Application certificate |
| `APPLE_INSTALLER_CERT_64` | Base64-encoded Developer ID Installer .p12 certificate |
| `APPLE_INSTALLER_PW` | Password for the Installer certificate |
| `APPLE_NOTARIZATION_PASSWORD` | App-specific password for notarization |
| `KEYCHAIN_PASSWORD` | Password for the temporary CI keychain |

The following are configured directly in the workflow (not secrets):
- **Signing Identity**: `Developer ID Application: OpenMS Inc. (C64UCGJ5PL)`
- **Installer Identity**: `Developer ID Installer: OpenMS Inc. (C64UCGJ5PL)`
- **Email**: `apple@openms.de`
- **Team ID**: `C64UCGJ5PL` (extracted from the identity string)

### Exporting Certificates

1. Open **Xcode → Settings → Accounts**
2. Select your Team and click **Manage Certificates**
3. Create or select **Developer ID Application** certificate
4. Right-click → **Export Certificate** → save as `.p12` with a password
5. Base64 encode for CI: `base64 -i certificate.p12 | pbcopy`
6. Repeat for **Developer ID Installer** certificate

### Finding Your Identity Strings and Team ID

Run this command to see your available signing identities:

```shell
security find-identity -v -p codesigning
```

Output example:
```
1) ABC123XYZ... "Developer ID Application: OpenMS Inc. (C64UCGJ5PL)"
2) DEF456ABC... "Developer ID Installer: OpenMS Inc. (C64UCGJ5PL)"
   2 valid identities found
```

- The full quoted string is your identity (e.g., `Developer ID Application: OpenMS Inc. (C64UCGJ5PL)`)
- The 10-character code in parentheses is your Team ID (e.g., `C64UCGJ5PL`)

### Creating an App-Specific Password

For notarization, create an app-specific password:

1. Go to https://appleid.apple.com
2. Sign in → Security → App-Specific Passwords → Generate
3. Use this password for `APPLE_NOTARIZATION_PASSWORD`

### Local Testing

To test code signing locally:

```shell
# Sign a binary
codesign --force --options runtime --timestamp -s "Developer ID Application: Your Name (TEAMID)" /path/to/binary

# Verify signature
codesign -dv --verbose=4 /path/to/binary

# Notarize a package (requires CODESIGNPW environment variable)
export CODESIGNPW="your-app-specific-password"
./cmake/MacOSX/notarize.sh package.pkg de.openms your@email.com CODESIGNPW .

# Verify notarization
spctl -a -vvv -t install package.pkg
```

### Troubleshooting

**"The specified item could not be found in the keychain"**
- Ensure the certificate is imported in the keychain
- Verify the identity string matches exactly (run `security find-identity -v`)

**Notarization fails with "Invalid"**
- Check if all binaries inside the package are signed with hardened runtime
- Get detailed log: `xcrun notarytool log <submission-id> --apple-id EMAIL --password PW --team-id TEAM log.json`

**"Developer cannot be verified" warning for users**
- The package was not notarized or stapling failed
- Users can override: System Settings → Privacy & Security → Allow

### Files in this Directory

- `notarize.sh` - Script for notarizing packages (.pkg, .dmg, .zip) using notarytool
- `sign_bins_and_libs.rb` - Ruby script to sign all binaries in a directory
- `fix_dependencies.rb` - Script to fix library paths for distribution
