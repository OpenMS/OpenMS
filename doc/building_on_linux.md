## Simplified Ubuntu Build Instructions (Recommended for Ubuntu ≥ 22.04)

These simplified instructions use only Ubuntu's official system packages and do not require building contrib libraries manually.

### Basic Installation (System Packages Only)

Include the Ubuntu "universe" repository and install dependencies:

`sudo add-apt-repository universe`

`sudo apt update`
`sudo apt-get install build-essential cmake` `autoconf patch libtool git automake`
`qt6-base-dev libqt6svg6-dev libeigen3-dev` `libboost-random-dev libboost-regex-dev`
`libboost-iostreams-dev libboost-date-time-dev`
 `libboost-math-dev libxerces-c-dev`
`zlib1g-dev libsvm-dev libbz2-dev` `coinor-libcoinmp-dev libhdf5-dev`


Clone the OpenMS repository:

`git clone https://github.com/OpenMS/OpenMS.git`


Configure and build OpenMS (recommended out-of-source build):

`mkdir OpenMS-build && cd OpenMS-build`
`cmake -DBOOST_USE_STATIC=OFF ../OpenMS`
`make -j4` # Adjust number of jobs according to CPU cores available


Set environment variables by adding to `~/.bashrc`:

`export LD_LIBRARY_PATH="/PATH/TO/OpenMS/lib:$LD_LIBRARY_PATH"`
`export PATH="$PATH:/PATH/TO/OpenMS/bin"`


Apply changes immediately:

`source ~/.bashrc`

### Advanced Option: Manual Contrib Build (Only if Necessary)

If your system lacks compatible package versions, you can manually build contrib libraries as follows:

Initialize contrib submodule and create build directory:

`cd ~/Development/OpenMS`
`git submodule update --init contrib`
`cd ..`
`mkdir contrib-build && cd contrib-build`
`cmake -DBUILD_TYPE=ALL -DNUMBER_OF_JOBS=4 ../OpenMS/contrib`


Then configure OpenMS using your built contrib libraries explicitly:

`cd ../OpenMS-build` # assuming you created this earlier or create again if needed
`cmake -DOPENMS_CONTRIB_LIBS="/PATH/TO/contrib-build" -DBOOST_USE_STATIC=OFF ../OpenMS`
`make -j4`



