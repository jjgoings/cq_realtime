# cq_realtime 

This is a python script to generate the absorption spectrum from a real time time-dependent SCF calculation in the output of [Chronus Quantum](https://github.com/liresearchgroup/chronusq_public).

## Use

You can run the script as-is on the provided output (water/Restricted Hartree-Fock/STO-3G), just:

```bash
$ python cq_realtime.py
```

The files are hard-coded into the main part of the script. You can always import the function `cqRealTime` in your own python scripts like so:
```python
from cq_realtime import cqRealTime
``` 

### Arguments
```python
def cqRealTime(real_time_file,dipole_direction,kick_strength,damp_const):
    '''
        CQ_RealTime.py: a post-processing script for computing the absorption spectrum of
         Real Time Time Dependent SCF jobs in Chronus Quantum

        Computes the energy range, w (in eV) and dipole strength function S(w) for
         a given real time TD-SCF run.

        real_time_file   ... type:string ; the RealTime_Dipole.csv file from a ChronusQ run
        dipole_direction ... type:char   ; which dipole moment contribution is computed (e.g. 'x','y', or 'z')
        kick_strength    ... type:float  ; in a.u., what was the applied field strength (e.g. 0.0001 au)
        damp_const       ... type:float  ; in a.u. of time, gives FWHM of 2/damp_const
    '''
```

Returns `w`, and `S`, which are the frequencies and absorption cross sections.

## Theory(-ish)

Three input files are included for water in a STO-3G basis. If you have ChronusQ installed, you can run these files.

To get the absorption spectrum you need to give the system three independent "perturbations," so we give the system a delta-function electric field "kick" at the beginning (e.g. turn on the field for one time step) in three different directions (x,y,z).

Then we measure the time-evolving dipole moment of the system along the polarization of that field. A Fourier transform of this time evolution is proportional to the dipole strength function, and also related to the absorption cross section (which we measure in lab in a UV-Vis experiment).






## Requirements

This code is written in Python 2.7 and requires [numpy](http://www.numpy.org/), [matplotlib](http://matplotlib.org/), and [SciPy](http://www.scipy.org/). On Ubuntu you can just: 

```bash
$ sudo apt-get install python-numpy python-scipy python-matplotlib 
```

Or on Fedora:
```bash
$ sudo yum install numpy scipy python-matplotlib 
```

For Mac OS X, if you use [Macports](http://www.macports.org/):
```bash
$ sudo port install py27-numpy py27-scipy py27-matplotlib 
```
