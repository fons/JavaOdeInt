#Summary

JavaOdeInt is a Java interface to a set of widely used Fortran libraries

# Building

## Dependencies

+ [maven](https://maven.apache.org/) of version 3.1 or later

+ [gfortran](https://gcc.gnu.org/wiki/GFortran) needs to be installed in order to build the source code

+ [bridj](https://github.com/nativelibs4java/BridJ) is used to generate the java interfaces. It's found in pom.xml in the jodeint directory.

+ [os-maven-plugin] (https://github.com/trustin/os-maven-plugin) The maven plugin used to determine the version of the operating system requires a maven version at least version 3.1

## Directory structure

The main directories
    
   + fodeint
   
     Each subdirectory corresponds to one of the Fortran packages [listed below](id:fortranpackages).  This contains the fortran code. Each directory name refers to the fortran source code package. Each package directory has two sub-directories :
     * src
       
       Contains the Fortran source code.
     * txt
      
      Contains the comment section of the Fortran source code describing the usage of the Fortran code as well as other auxiliary code or driver scripts where appropriate.
        
   
   + codeint
   
     This directory contains the C code interface. Each directory is the name of the fortran package it interfaces to preceded by a c. 
     E.g. the directory named *codepack* contains the c interface to *odepack* Fortran package.

## Mac OSX
  
  
  
## Linux


## Prebuilt



#[Fortran libraries covered](id:fortranpackages)


* [odepack](https://computation.llnl.gov/casc/odepack/odepack_home.html)

   The following routines are available :

    * lsoda
    * lsodar
    + lsode
    + lsodes


+ [dvode](https://computation.llnl.gov/casc/odepack/odepack_home.html)

+ [zvode](https://computation.llnl.gov/casc/odepack/odepack_home.html)

+ [quadpack](https://people.sc.fsu.edu/~jburkardt/f_src/quadpack/quadpack.html)

+ [radau5](http://www.unige.ch/~hairer/software.html)

+ [dopri5](http://www.unige.ch/~hairer/software.html)

+ [dop853](http://www.unige.ch/~hairer/software.html)

+ [gnicodes](http://www.unige.ch/~hairer/software.html)

+ [rkf45](https://people.sc.fsu.edu/~jburkardt/f77_src/rkf45/rkf45.html)


#Copyright

