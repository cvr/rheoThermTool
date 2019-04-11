#!/bin/bash
# vim: set fileencoding=utf-8 fileformat=unix :
# -*- coding: utf-8 -*-

cat << EOF
»»»»»»»»»»»»»»»»  rheoThermTool  ««««««««««««««««
»»»  Installation script for OpenFOAM-6  «««


This script will install rheoThermTool on a system with OpenFOAM-6

The information that the user should supply is:
WM_PROJECT_DIR        path for OpenFOAM-6 installation
RHEOTHERM_ROOT        path where rheoThermTool will be compiled
EIGEN_ROOT            path where to install Eigen
EOF

echo -e "What is the path for WM_PROJECT_DIR?"
read -e userInputProjectDir

if [ $? != 0 ] ; then
  echo "Error occurred while reading WM_PROJECT_DIR... aborting"
  exit 1
fi

eval userInputProjectDir=$userInputProjectDir

if [ ! -d $userInputProjectDir ] ; then
  echo "WM_PROJECT_DIR is not a directory... aborting"
  exit 1
elif [ ! -f $userInputProjectDir/etc/bashrc ] ; then
  echo "Failed to find WM_PROJECT_DIR/etc/bashrc... aborting"
  exit 1
fi

## Sourcing OpenFOAM paths
. $userInputProjectDir/etc/bashrc

if [ "$(realpath $userInputProjectDir)" != "$(realpath $WM_PROJECT_DIR)" ] ; then
  echo "Error: user input does not coincide with true WM_PROJECT_DIR"
  echo "       userInputProjectDir=${userInputProjectDir}"
  echo "            WM_PROJECT_DIR=${WM_PROJECT_DIR}"
  echo "       ... aborting"
  exit 1
fi


echo -e "In what path rheoThermTool source code will be compiled?"
echo -e "   (... to install at WM_THIRD_PARTY_DIR press [ENTER])"
read -e RHEOTHERM_ROOT

if [ $? != 0 ] ; then
  echo "Error occurred while reading RHEOTHERM_ROOT... aborting"
  exit 1
fi

eval RHEOTHERM_ROOT=$RHEOTHERM_ROOT

[ -z "${RHEOTHERM_ROOT}" ] && RHEOTHERM_ROOT=$WM_THIRD_PARTY_DIR
RHEOTHERM_ROOT=$(realpath $(readlink -f "${RHEOTHERM_ROOT}"))
echo "using RHEOTHERM_ROOT=${RHEOTHERM_ROOT}"

if [ ! -d $RHEOTHERM_ROOT ] ; then
  echo "RHEOTHERM_ROOT is not a directory... aborting"
  exit 1
elif [ ! -w $RHEOTHERM_ROOT ] ; then
  echo "Error: user has no write permissions on RHEOTHERM_ROOT... aborting"
  exit 1
fi


echo -e "In what path will Eigen be installed?"
echo -e "   (... to install at WM_THIRD_PARTY_DIR press [ENTER])"
read -e EIGEN_ROOT

if [ $? != 0 ] ; then
  echo "Error occurred while reading EIGEN_ROOT... aborting"
  exit 1
fi

eval EIGEN_ROOT=$EIGEN_ROOT

[ -z "${EIGEN_ROOT}" ] && EIGEN_ROOT=$WM_THIRD_PARTY_DIR
EIGEN_ROOT=$(realpath $(readlink -f "${EIGEN_ROOT}"))
echo "using EIGEN_ROOT=${EIGEN_ROOT}"

if [ ! -d $EIGEN_ROOT ] ; then
  echo "EIGEN_ROOT is not a directory... aborting"
  exit 1
elif [ ! -w $EIGEN_ROOT ] ; then
  echo "Error: user has no write permissions on EIGEN_ROOT... aborting"
  exit 1
fi



cd $RHEOTHERM_ROOT
git clone https://github.com/cvr/rheoThermTool rheoThermTool
cd rheoThermTool/
git log --decorate --graph > README_gitlog.md
VersionDate=$(git log -1 --format="%at" | xargs -I{} date -d @{} +%Y%m%d)
rm -rf .git .gitignore of+1806
cd ../
mv rheoThermTool "rheoThermTool-${VersionDate}"
ln -sf "rheoThermTool-${VersionDate}" rheoThermTool


## Install Eigen
cd $RHEOTHERM_ROOT/rheoThermTool
sed -i 's/WM_PROJECT_USER_DIR\/ThirdParty/WM_PROJECT_USER_DIR/g' downloadEigen
#sed -i 's/WM_PROJECT_USER_DIR/WM_THIRD_PARTY_DIR/g' downloadEigen
sed -i 's/WM_PROJECT_USER_DIR/EIGEN_ROOT/g' downloadEigen
export EIGEN_ROOT
./downloadEigen
EIGEN_RHEO=$(find "${EIGEN_ROOT}/" -maxdepth 1 -type d -iname "eigen*" | head -n 1)
export EIGEN_RHEO


## Create an environment variables script
cat <<EOF > $RHEOTHERM_ROOT/envars_rheoThermTool.sh
#!/bin/bash
#

## rheoThermTool
function RTT {
  export RHEOTHERM_DIR=$RHEOTHERM_ROOT/rheoThermTool
  export EIGEN_RHEO=$EIGEN_RHEO
  export RHEOTHERM_APPBIN=\$RHEOTHERM_DIR/platforms/\$WM_OPTIONS/bin
  export RHEOTHERM_LIBBIN=\$RHEOTHERM_DIR/platforms/\$WM_OPTIONS/lib
  export RHEOTHERM_TUTORIALS=\$RHEOTHERM_DIR/of60/tutorials
  export LD_LIBRARY_PATH=\$RHEOTHERM_LIBBIN:\$LD_LIBRARY_PATH
  export PATH=\$RHEOTHERM_APPBIN:\$PATH
}
RTT

echo -e "
rheoThermTool installation
~~~~~~~~~~~~~~~~~~~~~~~~~~
RHEOTHERM_DIR       = \$RHEOTHERM_DIR
RHEOTHERM_APPBIN    = \$RHEOTHERM_APPBIN
RHEOTHERM_TUTORIALS = \$RHEOTHERM_TUTORIALS
"

EOF
chmod +x $RHEOTHERM_ROOT/envars_rheoThermTool.sh


echo "Installing rheoThermTool..."

. $RHEOTHERM_ROOT/envars_rheoThermTool.sh

cd $RHEOTHERM_DIR/of60/src

env \
    WM_PROJECT_USER_DIR="${RHEOTHERM_DIR}" \
    FOAM_USER_APPBIN="${RHEOTHERM_APPBIN}" \
    FOAM_USER_LIBBIN="${RHEOTHERM_LIBBIN}" \
    ./Allwmake 2>&1 | tee -a log.Allwmake1


echo "Checking for compilation errors"
grep --color=auto -i error $RHEOTHERM_DIR/of60/src/log.Allwmake1


echo "Test rheoThermTool installation"
$RHEOTHERM_APPBIN/rheoThermFoam
rheoThermFoam


cat << EOF
To verify the installation, the user should run one of
the tutorials, e.g. the 2D channel.  For this, run:


source $userInputProjectDir/etc/bashrc
source $RHEOTHERM_ROOT/envars_rheoThermTool.sh

mkdir -p \$FOAM_RUN

run

cp -r \$RHEOTHERM_TUTORIALS/channel rheoThermFoam-channel

cd rheoThermFoam-channel

cp -r 0.orig 0

blockMesh

rheoThermFoam

paraFoam
EOF

