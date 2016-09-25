#Summary

JavaOdeInt is a Java interface to a set of widely used fortran libraries

# Building

## Dependencies

+ [maven](https://maven.apache.org/) of version 3.1 or later

+ [gfortran](https://gcc.gnu.org/wiki/GFortran) needs to be installed in order to build the source code

+ [bridj](https://github.com/nativelibs4java/BridJ) is used to generate the java interfaces. It's found in pom.xml in the jodeint directory.

+ [os-maven-plugin] (https://github.com/trustin/os-maven-plugin) The maven plugin used to determine the version of the operating system requires a maven version at least version 3.1


## Mac OSX
  
  
  
## Linux

## Prebuilt



#Fortran libraries covered


* [odepack](https://computation.llnl.gov/casc/odepack/odepack_home.html)

   The following routines are available :

    * lsoda
    * lsodar
    + lsode
    + lsodes


+ dvode

+ zvode

+ quadpack

+ radau5

+ dopri5

+ dop854

+ gnicodes

+ rkf45


#Copyright

