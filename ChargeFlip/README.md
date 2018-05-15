# Usage

+ ### Check out from GitLab

+ By ssh 

+ git clone ssh://git@gitlab.cern.ch:7999/hku/DataDriven.git

+ Or by https

+ git clone https://gitlab.cern.ch/hku/DataDriven.git

+ ### Create ROOT-based environment

+ cd DataDriven

+ setupATLAS

+ rcSetup Base 2.4.31

+ rc find_packages

+ rc clean

+ ### Complie and create some symbol links

+ rc compile

+ cd ChargeFlip/run

+ ln -s ../../RootCoreBin/bin/x86_64-slc6-gcc49-opt/charge_flip charge_flip

+ ln -s ../../RootCoreBin/bin/x86_64-slc6-gcc49-opt/re_weight re_weight

+ ### Run util/charge_flip.cxx

+ ./run_data.sh (for data)

+ ./run_mc.sh   (for MC)

+ ### Run util/re_weight.cxx

+ ./re_weight_data.sh (for data)

+ ./re_weight_mc.sh   (for MC)

+ ### Perform fit
