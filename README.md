[![Build Status](https://api.travis-ci.org/predictionmachines/Filzbach.svg)](https://travis-ci.org/predictionmachines/Filzbach)
[![Build status](https://ci.appveyor.com/api/projects/status/b5w7f2wcvgtj24bw?svg=true)](https://ci.appveyor.com/project/vassilyl/filzbach)


#Filzbach

_Fit complex models to heterogeneous data: Bayesian and likelihood analysis made easy_

Filzbach is a flexible, fast, robust, parameter estimation engine
that allows you to parameterize arbitrary, non-linear models,
of the kind that are necessary in biological sciences, against multiple, heterogeneous data sets.
Filzbach allows for Bayesian parameter estimation, maximum likelihood analysis,
priors, latents, hierarchies, error propagation and model selection, from just a few lines of code.

Home web site: http://research.microsoft.com/filzbach


[Filzbach User Guide](http://research.microsoft.com/en-us/um/cambridge/groups/science/tools/filzbach/Filzbach%20User%20Gude%20v.1.1.pdf)

#Build Filzbach library

##Windows

To build static library `filzbach.lib` in "Debug" configuration from command prompt:

```
cd filzbach
msbuild
```

To build a dynamically linked library `filzbach.dll` in "Release" configuration:

```
cd filzbach
msbuild /p:Configuration=Release_DLL
```

To build either a static or a dynamic library with support for OpenMP set the CL environment variable which affects Visual C++ compiler options, e.g.

```
set CL=/openmp
cd filzbach
msbuild /p:Configuration=Release_DLL
```
 

##Unix

Build a library `libfilzbach.a` from a shell prompt: 

```
cd filzbach
make
```

#Run examples

The repository contains several sample programs that illustrate how to apply Filzbach to different kinds of models.

Id | Name | Desciption
---|------|-------------
1 | normal | How to learn parameters of a normal distribution
2 | normal_priors | How to set up priors
3 | mixednormal | Learn the location and width of several peaks 
4 | poisson | Learn parameters of a Poisson distribution
5 | poisson_multispp| Learn parameters of several Poisson distributions at once with hierarchical modelling
7 | lr | Linear regression
8 | logistic_regression | Logistic regression
9 | stepfunction | Step function
11| plantgrowth | Infer parameters of a simple plant growth model from observations
12| seedrain | Seed 
13| SDM | Species distribution modelling
16| ProteinInteraction | Simple model of a protein interaction network
17| TreeMort | A model of tree mortality

On Windows, the easiest is to start [Microsoft Visual Studio](http://microsoft.com/visualstudio), 
open filzbach.sln, uncomment one of the lines in examples/examples.h header file and run the examples project.

You can also compile and run from command prompt. You then first set set up the CL environment variable with an example Id. 
E.g., to run species distribution modelling sample (Id=13) do the following:

```
set CL=/DMODEL#13
msbuild /t:Clean
msbuild
Debug\examples
```

On Unix systems from the `examples` directory use `make` to build an individual example using its name or build all the examples with one command:

```
cd examples
make examples
```
