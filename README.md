Raremetal2: A tool for rare variants meta-analysis
	(c) 2012-2017 Sai Chen, Shuang Feng, Dajiang Liu, Gonçalo Abecasis

For more info, please see wiki page
[https://genome.sph.umich.edu/wiki/RAREMETAL](https://genome.sph.umich.edu/wiki/RAREMETAL)

Raremetal2 is an updated version of the meta-analysis software Raremetal
It is 30% faster than the original version and takes less memory
Compatible with old Raremetalworker or rvtest generated summary files

This is still a beta version. Single-variant meta-analysis looks good on T2D and MGI data,  but needs more 
optimizations on group test.

New features:
1. new method for unbalanced studies. To toggle that, please specify --useExact
2. multi-allelic site meta-analysis enabled. However, if you have Raremetalworker or rvtest covariance files from 
  previous versions, multi-allelic site meta-analysis won't work well in group test because these old covariance files 
  do not record alleles for each site.
3. missing data in conditional analysis


## Installation
### Requirements
Raremetal2 requires `libRMath` and `zlib` to be installed on your system.

### Build
From the project root, run `make` to build binaries (will be deposited in `bin/`).

Run `make clean` to remove build assets from a previous attempt. 

Currently, only Linux operating systems are officially supported.

