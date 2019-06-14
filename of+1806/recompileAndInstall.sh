#!/bin/bash
# vim: set fileencoding=utf-8 fileformat=unix :
# -*- coding: utf-8 -*-

## Source OpenFOAM and rheoThermTool environment variables
. /opt/OpenFOAM/OpenFOAM-v1806/etc/bashrc

cd $RHEOTHERM_DIR/of+1806

while true; do
	read -p "Update code with that one in the repository [y-yes, n-no]? " yn
	case $yn in
		[Yy]*)
			git clone https://github.com/cvr/rheoThermTool.git
			if [ $exitCondition -ne 0] ; then
				echo "Error: failed to clone rheoThermTool, aborting"
				rm -rf rheoThermTool
				exit 1
			fi
			rsync -Pa --delete rheoThermTool/of+1806/tutorials ./
			rsync -Pa --delete rheoThermTool/of+1806/src ./
			rm -rf rheoThermTool
			break
			;;
		[Nn]*)
			break
			;;
		*)
			echo "Please, answer yes or no"
			;;
	esac
done

## Recompile
cd $RHEOTHERM_DIR/of+1806/src

env \
	WM_PROJECT_USER_DIR="${RHEOTHERM_DIR}" \
	FOAM_USER_APPBIN="${RHEOTHERM_APPBIN}" \
	FOAM_USER_LIBBIN="${RHEOTHERM_LIBBIN}" \
	./Allwmake 2>&1 | tee -a log.Allwmake

