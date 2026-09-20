# Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: BSD-3-Clause
# $Maintainer: Timo Sachsenberg $

# PCH caching needs the MSVC invalidation/path fixes in ccache 4.14. Use the
# upstream binary so the opt-in CI mode does not depend on Chocolatey lag.
$ErrorActionPreference = 'Stop'
$destination = Join-Path $env:RUNNER_TEMP 'openms-pch-ccache'
$archive = Join-Path $env:RUNNER_TEMP 'ccache-4.14-windows-x86_64.zip'
Invoke-WebRequest -Uri 'https://github.com/ccache/ccache/releases/download/v4.14/ccache-4.14-windows-x86_64.zip' -OutFile $archive
if ((Get-FileHash $archive -Algorithm SHA256).Hash -ne '2568347a697e103ca1b073981c704ad76fb2507d066c38dba038dd73399d968f') {
    throw 'ccache archive checksum mismatch'
}
Expand-Archive -Path $archive -DestinationPath $destination -Force
$binaries = @(Get-ChildItem $destination -Filter ccache.exe -File -Recurse)
if ($binaries.Count -ne 1) {
    throw 'Expected exactly one ccache executable'
}
$binaries[0].DirectoryName | Out-File -FilePath $env:GITHUB_PATH -Encoding utf8 -Append
& $binaries[0].FullName --version
if ($LASTEXITCODE -ne 0) {
    throw 'ccache executable failed'
}
