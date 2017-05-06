# CharFun: The Characteristic Functions Toolbox
MATLAB repository of characteristic functions and tools for their combinations and numerical inversion.

For the R version of the toolbox see the CharFun package development available at

- https://github.com/Simkova/CharFun


For current status of the MATLAB toolbox (not an identical clone) see the CharFunTool development available at

- https://github.com/witkovsky/CharFunTool

About
=====

The Characteristic Functions Toolbox (CharFun) consists of a set of algorithms for evaluating selected characteristic funcions
and algorithms for numerical inversion of the (combined and/or compound) characteristic functions, used to evaluate the probability density function (PDF) and the cumulative distribution function (CDF).
                                                                              
The toolbox includes inversion algorithm, including those based on simple trapezoidal rule for computing the integrals defined by the Gil-Pelaez formulae, and/or by using the FFT algorithm for computing the Fourier transform integrals.
                                                                       
Installation and requirements
=============================

CharFun was developed with R version 3.3.1 (2016-06-21).

To install, you can either clone the directory with Git or download a .zip file or install package from R CRAN-repository.

## Option 1: Download .zip file

Download a .zip of CharFun from

- https://github.com/Simkova/CharFun/releases

After unzipping, you will need open CharFun.Rproj.

## Option 2: Clone with Git

To clone the CharFun repository, first navigate in a terminal to where you want the repository cloned, then type
```
git clone https://github.com/simkova/CharFun.git
```
and you will need open CharFun.Rproj.

## Option 3: Install package

You can download package from 

- https://github.com/Simkova/CharFun/releases

in installe packages in R-studio you chose Package Archive File (.tar.gz).


Getting started
===============

We recommend taking a look at the Examples collection. 

To get a taste of what computing with CharFun is like, type
```
   cf <- function(t) exp(-t^2/2)  # the standard normal characteristic function (CF)
   result <- cf2DistGP(cf)   # Invert the CF to get the CDF and PDF   
```


License
=======

See `LICENSE.txt` for CharFunTool licensing information.
