# What VMStools can do for you #

## Introduction ##

Vmstools is an open-source software package build in R. It is specifically developed to process, analyse and visualize logbooks and VMS data.

### Details ###
Some of its functionality:
  * Analyses use standardised data formats for logbook (EFLALO) and VMS (TACSAT)
  * Handles erroneous data point detection and removal
  * Métier identification through the use of clustering techniques
  * Link logbook and VMS data together in order to distinguish fishing from other activities
  * Provide high-resolution maps of both fishing effort and -landings, interpolate vessel tracks
  * Calculate indicators of fishing impact as listed under the Data Collection Framework at different spatio-temporal scales.
  * Prepare data for FishFrame submission
  * Link other spatial data to VMS and logbook data

## System Requirements ##
You will need the following software to use the VMStools software and repository.

  * R (version 2.14.x or 2.15.x) for your system. Available [here](http://www.r-project.org).
  * The VMStools package itself, available for download [here](http://code.google.com/p/vmstools/downloads/)
  * For very enthusiastic users / developers: An SVN (subversion) client. For Windows users, we recommend TortoiseSVN. Linux users will typically have access to such a client through their package distribution system.
  * Some way to write and edit R code. For Windows, try RStudio or Tinn-R. For Linux, you almost certaintly already have Vi

## Getting the Code for enthusiastic users / developers ##
Once you've got all the necessary bits and pieces, the next step is to fetch a "working copy" of the VMStools repository from the server. The SVN client that you have installed takes care of this. The general process is client specific, but generally involves two steps
  1. Create a directory for your working copy (e.g. Desktop/VMStools) and
  1. "Checking out" the repository.

We'll assume you can handle 1 and 2. yourself, but number 3. is a bit trickier, and is of course dependent on your SVN client. There are a number of good resources available e.g. for windows

Firstly there are two different ways to checkout the code, depending on whether you want to commit the changes back to the repository or not. If you checkout the repository from http://vmstools.googlecode.com/svn/trunk/ it will not be possible to commit your changes back - however, the benefit is that anyone can get access to the repository, without the need for a google account.
Alternatively, you can check out the code from https://vmstools.googlecode.com/svn/trunk/ (note the https, rather than http), in which case you will be able to commit your changes back. The disadvantage is that this requires a google account (e.g. gmail.com).
Finally, there is a wee trap for beginners when logging in for the first time - make sure that you use the password for the project that googlecode supplies you - not your login password to google itself (slightly counter-intuitive, I know)!

## Installing the VMStools package ##
Download the VMStools zip or rar file from the repository under the download tab or build it yourself from the repository. Please note that the VMStools package depends on a number of other R packages. They are listed in the DESCRIPTION file in /trunk of the repository.
The packages VMStools depends on can easily be installed through:
```
vmstoolsPackages <-  c("cluster","data.table","doBy","maps","mapdata","maptools","PBSmapping","sp")

for(i in vmstoolsPackages)       try(install.packages(pkgs=i,repos=getOption("repos")))```

Once you've installed all secondary packages install VMStools using ```
install.packages("yourPath/vmstools.zip")```