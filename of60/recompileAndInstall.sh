#!/bin/bash
# vim: set fileencoding=utf-8 fileformat=unix :
# -*- coding: utf-8 -*-

## Source OpenFOAM and rheoThermTool environment variables
. /opt/OpenFOAM/OpenFOAM-6/etc/bashrc

## Recompile
cd $RHEOTHERM_DIR/of60/src

env \
	WM_PROJECT_USER_DIR="${RHEOTHERM_DIR}" \
	FOAM_USER_APPBIN="${RHEOTHERM_APPBIN}" \
	FOAM_USER_LIBBIN="${RHEOTHERM_LIBBIN}" \
	./Allwmake 2>&1 | tee -a log.Allwmake2

