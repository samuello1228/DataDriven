# Introduction

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

### Usage of ChargeFlip packages

+ [see this page](ChargeFlip/README.md)
