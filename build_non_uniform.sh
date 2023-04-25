REPO_ROOT="$(dirname -- "${BASH_SOURCE[0]}")"            # relative
REPO_ROOT="$(cd -- "$REPO_ROOT" && pwd)"    # absolutized and normalized
if [[ -z "$REPO_ROOT" ]] ; then
  # error; for some reason, the path is not accessible
  # to the script (e.g. permissions re-evaled after suid)
  exit 1  # fail
fi

ALPSCore_DIR="$REPO_ROOT"/alpscore
if ! [ -d "$ALPSCore_DIR" ]; then
    mkdir "$ALPSCore_DIR"
fi

if [ -d "$ALPSCore_DIR"/ALPSCore ]; then
    # skip
    echo "ALPSCore already exists"
else
    cd "$ALPSCore_DIR" && git clone https://github.com/ALPSCore/ALPSCore.git

    cd "$ALPSCore_DIR"/ALPSCore
    mkdir build
    cd build

    # if eigen headers are given, use them
    if ! [ -z "$1" ]; then
        cmake ../ -DCMAKE_INSTALL_PREFIX="$ALPSCore_DIR"/install -DEIGEN3_INCLUDE_DIR="$1"
    else
        cmake ../ -DALPS_INSTALL_EIGEN=yes -DCMAKE_INSTALL_PREFIX="$ALPSCore_DIR"/install -DEIGEN3_INCLUDE_DIR="$1"
    fi

    make -j 5
    make test
    make install
fi



BUILD_DIR="$REPO_ROOT"/build_non_uniform
if [ -d "$BUILD_DIR" ];
then
    rm -rf "$BUILD_DIR"
fi

mkdir "$BUILD_DIR"
cd "$BUILD_DIR" && cmake "$REPO_ROOT"/src -DCMAKE_PREFIX_PATH="$REPO_ROOT"/alpscore/install -DUNIFORM=False -DLATTICE=square && make -j5

