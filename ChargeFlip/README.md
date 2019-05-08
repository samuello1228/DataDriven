# Usage

### Check out from GitLab

+ By ssh 
	+ git clone ssh://git@gitlab.cern.ch:7999/hku/DataDriven.git
+ By https
	+ git clone https://gitlab.cern.ch/hku/DataDriven.git

### Create ROOT-based environment

+ cd DataDriven
+ setupATLAS
+ rcSetup Base 2.4.31
+ rc find_packages
+ rc clean

### Complie and create some symbol links

+ rc compile
+ cd $ROOTCOREBIN/../ChargeFlip/run
+ ln -s $ROOTCOREBIN/bin/x86_64-slc6-gcc49-opt/charge_flip charge_flip
+ ln -s $ROOTCOREBIN/bin/x86_64-slc6-gcc49-opt/re_weight re_weight

### Run util/charge_flip.cxx

+ cd $ROOTCOREBIN/../ChargeFlip/run
+ ./run_data.sh (for data)
+ ./run_mc.sh   (for MC)

### Perform fit

+ cd $ROOTCOREBIN/../ChargeFlip/scripts/Fit
+ root -l -b -q fit.c+

### Extract MC Truth

+ cd $ROOTCOREBIN/../ChargeFlip/scripts/FitPlots
+ root -l -b -q get_truth.cxx+

### Draw fitting results 

+ cd $ROOTCOREBIN/../ChargeFlip/scripts/FitPlots
+ root -l -b -q draw.c+

### Encapsulate the fitted results into ROOT files

+ cd $ROOTCOREBIN/../ChargeFlip/scripts/GenNTuples
+ root -l -b -q gen_signal.c+

### Run util/re_weight.cxx

+ cd $ROOTCOREBIN/../ChargeFlip/run
+ ./re_weight_data.sh (for data)
+ ./re_weight_mc.sh   (for MC)

### Draw plots of re-weighted control samples for validation

+ cd $ROOTCOREBIN/../ChargeFlip/scripts/ReweightPlots
+ root -l -b -q draw.c+

