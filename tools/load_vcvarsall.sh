# This is a script for creating a "Visual Studio Command Line" in (git) bash on Windows
# It basically invokes the corresponding vcvarsall.bat (that is also called when opening a VS Command Line)
# and copies all the environment variables from it.

# You can source it or add this function to your ~/.bashrc for example

# call load_msenv [year] from a shell to load the x64 enviroment for the corresponding VS year.
# it will cache the variables in $HOME/.msenv${msversion_year}_bash. Delete it if you want to reload.
function load_msenv() {
    local msversion_year=2022
    if [ $# -gt 0 ]; then
        msversion_year=$1
    fi

    case $msversion_year in
        2017)
            msversion_prod=15
            ;;
        2019)
            msversion_prod=16
            ;;
        2022)
            msversion_prod=17
            ;;
        *)
            >&2 printf "Invalid version year. Supported are 2017, 2019 and 2022."
            return 1
    esac

    msversion_prod_p1=$(($msversion_prod+1))

    local VSWHERE='C:\Program Files (x86)\Microsoft Visual Studio\Installer\vswhere.exe'
    installPath=$("$VSWHERE" -products '*' -version "[$msversion_prod,$msversion_prod_p1)" -property installationPath)

    # TODO check if exists, print error otherwise
    local vcvarsall="${installPath}\\VC\\Auxiliary\\Build\\vcvarsall.bat"

    local OLDPATH=$PATH
    local msenv="$HOME/.msenv_${msversion_year}_bash"
    if [ ! -f "$msenv" ]; then
        local msenvbatch="__print_ms_env.bat"
        echo "@echo off" > "$msenvbatch"
        echo "call \"${vcvarsall}\" x64" >> "$msenvbatch"
        echo "set" >> "$msenvbatch"
        cmd "/C $msenvbatch" > "$msenv.tmp"
        rm -f "$msenvbatch"
        
        grep -E '^PATH=' "$msenv.tmp" | \
            sed \
                -e 's/\(.*\)=\(.*\)/export \1="\2"/g' \
                -e 's/\([a-zA-Z]\):[\\\/]/\/\1\//g' \
                -e 's/\\/\//g' \
                -e 's/;\//:\//g' \
                > "$msenv"

        # Don't mess with CL compilation env
        grep -E '^(INCLUDE|LIB|LIBPATH)=' "$msenv.tmp" | \
            sed \
                -e 's/\(.*\)=\(.*\)/export \1="\2"/g' \
                >> "$msenv"

        rm "$msenv.tmp"
    fi

    source "$msenv"
    export PATH="$PATH:$OLDPATH"
}

export -f load_msenv