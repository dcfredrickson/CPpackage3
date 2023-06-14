# CPpackage

A software tool for performing DFT-Chemical Pressure analysis using abinit.

## Installation

Provided you have an SSH key to clone the repo, just run these commands:
```
git clone git@github.com:dcfredrickson/CPpackage
cd CPpackage
make
cd bin
```
This wil clone the repo into your current directory, change to it, run the 
makefile (which automatically compiles everything), and changes directories
again to get you to the binaries. You can then copy the binaries to any
directory in your `$PATH` (`/usr/local/bin`, for instance).

You can run `git pull` in the repo to get any future updates to CPpackage.

## Dependencies

### Compilation

The GNU Scientific Library and libxc are required for compiling CPpackage.

### Running calculations

Both CPpackage and the nonlocal programs take outputs from abinit version
7.10.5. The method requires the use of HGH (Hartwigsen-Goedecker-Hutter) 
norm-conserving pseudopotentials, which are available
[here.](https://www.abinit.org/psp-tables) The package currently supports the
LDA and PBE exchange-correlation functionals.

Atomic core profiles required as part of the Hirshfeld-inspired scheme may be
generated with 
[APE (Atomic Pseudopotentials Engine).](https://www.tddft.org/programs/APE/)

## Citations and more information

If you use CPpackage in your research, you can cite the following papers:

1. Lu, E.; Van Buskirk, J. S.; Cheng, J.; Fredrickson, D. C. Tutorial on 
Chemical Pressure Analysis: How Atomic Packing Drives Laves/Zintl Intergrowth in
K3Au5Tl. *Crystals* **2021**, *11* (8), 906. 
https://doi.org/10.3390/cryst11080906.

2. Berns, V. M.; Engelkemier, J.; Guo, Y.; Kilduff, B. J.; Fredrickson, D. C.
Progress in Visualizing Atomic Size Effects with DFT-Chemical Pressure
Analysis: From Isolated Atoms to Trends in AB5 Intermetallics. *J. Chem. Theory
Comput.* **2014**, *10* (8), 3380–3392. https://doi.org/10.1021/ct500246b.

3. Engelkemier, J.; Berns, V. M.; Fredrickson, D. C. First-Principles 
Elucidation of Atomic Size Effects Using DFT-Chemical Pressure Analysis: 
Origins of Ca36Sn23’s Long-Period Superstructure. *J. Chem. Theory Comput.* 
**2013**, *9* (7), 3170–3180. https://doi.org/10.1021/ct400274f.

4. Fredrickson, D. C. DFT-Chemical Pressure Analysis: Visualizing the Role of
Atomic Size in Shaping the Structures of Inorganic Materials. *J. Am. Chem. 
Soc.* **2012**, *134* (13), 5991–5999. https://doi.org/10.1021/ja300685j.
