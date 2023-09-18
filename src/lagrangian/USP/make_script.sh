#!/bin/sh

cd ~/OpenFOAM/nikovasi-v2112/src/lagrangian
./Allwmake

cd ~/OpenFOAM/nikovasi-v2112/applications
./Allwmake
