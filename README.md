#Summary

There is no need to reinvent the wheel when it comes to ODE solvers. 

JavaOdeInt is a Java interface to a set of widely respected Fortran ODE packages. 
These are the same packages powering Python's numpy library or  R's deSolve. 

JavaOdeInt 's aim is to provide the basics from which a fully fleshed out ode solver can be constructed.

Between the Java interface and the Fortran packages sits a very thin C-layer. 
The purpose of C interface is two-fold: To provide [bridj](https://github.com/nativelibs4java/BridJ) with a header file from which it can generate the Java classes and to provide a very simple interface to the Fortran routines.

This simple C interface only requires a callback function, dimensionality, the initial conditions and the time interval, plus any additional parameters. The results are written to a user provided stack.  

The C function executes a loop over the fortran function and writes the results to the user supplied  stack. It may also change some of the integration parameters if it detects that a larger number os steps are needed.




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

