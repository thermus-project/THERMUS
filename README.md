# THERMUS
A Thermal Model Package for ROOT v4

## Contacts

Thermus.Project@cern.ch
Thermus.Support@cern.ch

## TL-DR

### Inplace install and run
- Ensure ROOT (https://root.cern.ch) is installed and configured in your current shell.
- Build and configure for local run:
    ```bash
    ./scripts/inplace_build.sh
    ```
- Starts a root session with Thermus libs:
    ```bash
    ./run_thermus
    ```
## Particle datasets

### Base datasets

The Base dataset is installed under $THERMUS/share/Thermus/particles/base directory. It's the one used in older versions of THERMUS, before 4.0.

It contains the following subsets:
- PPB2002 
- PPB2014_CBHN
- PPB2014_CBHN_fixed_saveQM14
- PPB2018
- PPB2021

PPB2002 will likely not work in 4.0. # FIXME : ensure that.

### 2024 dataset

The 4.0 version comes with a new dataset in $THERMUS/share/Thermus/particles/2024 directory.

Its currents subset are:
- PPB2024

### Selecting a dataset and subset

To select one of those set and subsets, instantiate the TTMParticleSet the following way, for example for the 2024 PPB2024 subset:
```C++
TTMParticleSet myset("2024","PPB2024");
```
See next section for more details about custom sets and particles cascaded selection.

### Custom sets and particles cascaded selection

Particles set will be fetched from, in order, reading the first found:
- $CUSTOM_THERMUS_PARTICLES_DIR/particles/\<set>/\<subset>.txt
- $PWD/particles/\<set>/\<subset>.txt
- $THERMUS/share/Thermus/particles/\<set>/\<subset>.txt
- $THERMUS/share/Thermus/particles/base/\<subset>.txt

From a particle set, particles data will be fetched from, in order, reading the first found:
- $CUSTOM_THERMUS_PARTICLES_DIR/myparticles/\<set>/\<particle>.txt
- $PWD/myparticles/\<set>/\<particle>.txt
- $THERMUS/share/Thermus/particles/\<set>/\<particle>.txt
- $THERMUS/share/Thermus/particles/base/\<particle>.txt

So, for example, for the 2024 PPB2024 subset, if $CUSTOM_THERMUS_PARTICLES_DIR is set to "/home/myuser/custom_thermus_dir" and thermus is ran from "/home/myuser/myproj", then, the file "PPB2024.txt" will be fetched from:
- /home/myuser/custom_thermus_dir/particles/2024/PPB2024.txt
If not existing, it will be fetched from:
- /home/myuser/myproj/particles/2024/PPB2024.txt
If not existing, the default path will be used:
- $THERMUS/share/Thermus/particles/2024/PPB2024.txt
(and that file exists).

When the particle file "Sigma0_decay.txt" will be requested, it will again be read from either:
- /home/myuser/custom_thermus_dir/particles/2024/Sigma0_decay.txt
If not found, it will be read from:
- /home/myuser/myproj/particles/2024/Sigma0_decay.txt
If not found, it will be read from:
- $THERMUS/share/Thermus/particles/2024/Sigma0_decay.txt
And then, base will be used:
- $THERMUS/share/Thermus/particles/base/Sigma0_decay.txt

Note: the failover to base allows not to duplicate all files in the more recent sets. Only modified files are needed.

## Contents

Directory THERMUS/ :

README      - file containing important information
LICENSE     - usage terms and conditions
functions   - functions needed for numerical calculations
include     - prototypes for constrains and models
main        - models and main classes
particles   - list of particles and decay properties
test        - exemples and basic tests

## Environment variables

The environment variables ROOTSYS and THERMUS should be set and exported
respectively with the proper <pathnames> (checked by cmakethermus.sh):
    $ROOTSYS should point to the directory of the root installation
    $THERMUS should point to the directory containing the THERMUS code

## Compiling THERMUS and running test programs

Configuration and compilation should be done using CMake: source cmakethermus.sh
(Please note that the provided makefiles are not recommended).
For using THERMUS, a ROOT version >= 6 compiled from source is necessary:
you can obtain the source via: git clone http://root.cern.ch/git/root.git root
please look at cmakeroot.sh as an example for properly building and installing root
Note that a rootlogon.C loading THERMUS Shared Object (.so) libraries is also provided.


