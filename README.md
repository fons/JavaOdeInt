#Summary

There is no need to reinvent the wheel when it comes to ODE solvers. 

JavaOdeInt provides Java bindings to a set of widely respected Fortran ODE packages. 
These are the same packages powering Python's numpy library or  R's deSolve. 

JavaOdeInt 's aim is to provide the basics from which a fully fleshed out ode solver can be constructed.

Between the Java interface and the Fortran packages sits a very thin C-layer. 
The purpose of C interface is two-fold: To provide [bridj](https://github.com/nativelibs4java/BridJ) with a header file from which it can generate the Java classes and to provide a very simple interface to the Fortran routines.

This simple C interface requires only the very basics :a callback function, the dimensionality of the ode system, the initial conditions and the time interval, plus any additional parameters. The results are written to a user provided stack.  

The C function executes a loop over the fortran function and writes the results to the user supplied  stack. It may also change some of the integration parameters if it detects that a larger number os steps are needed.

Using the full Fortran interface is also possible. Each Fortran function has an extensive preamble describing its use. This was extracted from the Fortran package and put in its txt subdirectory.

Regardless of the interface you need to use [bridj](https://github.com/nativelibs4java/BridJ) for memory allocation and pointer management.

Both the C interface and the Java interface come with a set of tests. The tests include examples of the use of the simple and full interface.

# Building

## Dependencies

+ [maven](https://maven.apache.org/) of version 3.1 or later

+ [gfortran](https://gcc.gnu.org/wiki/GFortran) needs to be installed in order to build the source code

+ [bridj](https://github.com/nativelibs4java/BridJ) is used to generate the java interfaces. It's found in pom.xml in the jodeint directory.

+ [os-maven-plugin] (https://github.com/trustin/os-maven-plugin) The maven plugin used to determine the version of the operating system requires a maven version at least version 3.1

## Directory structure

The main directories
    
   + __fodeint__ : The Fortran modules
   
     Each subdirectory corresponds to one of the Fortran packages [listed below](id:fortranpackages) and contains the Fortran source code. Each package directory has two subdirectories. The __src__ subdirectory contains the Fortran source code and Makefile.
The  __txt__ subdirectory contains the preamble of the Fortran source code. It also contains any available  auxiliary code or driver scripts.     
     
 
   + __codeint__ : The glue
   
     This directory contains the C code interface. The name of each subdirectory is the name of the fortran package it interfaces to prefixed by a c.
      So the directory named *codepack* contains the c interface to *odepack* Fortran package. 
      
      Each sub directory  has two more directories :  __src__ and __include__. __src__ contains the source code, internal headers and makefile.  __include__  contains the header used by [bridj](https://github.com/nativelibs4java/BridJ) to generate the Java interface classes. 

   + __jodeint__ : The Java interface
    
   This directory contains the maven pom.xml and the Java code for the tests. The [bridj](https://github.com/nativelibs4java/BridJ) generator will create the code for the Java bindings on the fly. These are found in the __src__ path and can be inspected. You are not expected to edit these files in any way.
   
 The __target__ directory contains files generated by the maven build process.

 
 
   + config
   
   Contains the main make files used by the Fortran and C code.
   
   + data
   
   Contains the reference data for the Java regression tests.
   
   + util 
   
   Contains various utilities.
 
   + target
   
   This will contain the Fortran and C object and library files. 
   
   For the Linux version you may need to set LD_LIBRARY_PATH to the target library path. 


## Mac OSX
  
  mvn clean package
  
## Linux

    export LD_LIBRARY_PATH=`pwd`/src/main/resources/lib/linux_x64/:$LD_LIBRARY_PATH
    mvn clean package

## Prebuilt

 TBA
 
## Build process

The main build script is __Build.sh__ which can be found in ./src/scripts in the jodeint folder.

The build script executes the make file for each of the Fortran and C code packages, copies the headers and libraries to locations maven knows about, and then uses maven to run the  [bridj](https://github.com/nativelibs4java/BridJ) generator on each of the C interface header files. It does this by creating a configuration script for each package and hard linking bridj's configuration script to each one of these library specific ones.



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

