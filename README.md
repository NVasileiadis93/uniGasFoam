# Prerequisites
* openFOAM-v2112 or later
* python3.8 or later

# Compiling hybridDCFoam
[//]: # (source openFOAM bashrc)
* source ~/OpenFOAM/OpenFOAM-**[VERSION]**/etc/bashrc

[//]: # (compile hybridDCFoam)
* ./Allwmake [nprocs]

[//]: # (create symbolic link)
* sudo ln -s $WM_PROJECT_USER_DIR/scripts/hybridDCFoam/dist/hybridDCFoam /usr/local/bin/hybridDCFoam
