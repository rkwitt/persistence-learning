# Persistence Learning


This code implements the **kernel(s) for persistence diagrams** proposed in the following
two publications. Please use the provided BibTeX entries when citing our work.

```
@inproceedings{Reininghaus14a,
    author    = {R.~Reininghaus and U.~Bauer and S.~Huber and R.~Kwitt},
    title     = {A Stable Multi-scale Kernel for Topological Machine Learning},
    booktitle = {CVPR},
    year      = {2015}}
```
```
@inproceedings{Kwitt15a,
    author    = {R.~Kwitt and S.~Huber and M.~Niethammer and W.~Lin and U.~Bauer},
    title     = {Statistical Topological Data Analysis - A Kernel Perspective},
    booktitle = {NIPS},
    year      = {2015}}
```

# Overview

**[Compilation](#compilation)**   
**[Examples](#examples)**

# Compilation

The core of the code is in ```diagram_distance.cpp``` and depends on
DIPHA which is included as a submodule. After you have checked out the
repository via

```bash
git clone https://github.com/rkwitt/persistence-learning.git
```

you can checkout the submodule(s) via

```bash
git update submodule --recursive
```

In case the submodules do not get checked-out properly (or nothing happens),
execute the heper script

```bash
./git-submodule-sync.rb
```

Once this has finished, change into the ```code/diagram_distance``` directory
and create a ```build``` directory, then use ```cmake``` or ```ccmake``` to
configure the build process, e.g.,

```bash
cd code/dipha-pss
mkdir build
cd build
cmake ...
make
```

If you want to build ```diagram_distance``` *without* support for the
Fast-Gauss-Transform, simply take the standard settings as they are. In
case you do *not* have MPI support on your machine, you can install, e.g.,
OpenMPI. On MacOS (using *homebrew*) this can simply be done via

```bash
brew install open-mpi
```

### Compiling support for the Fast-Gauss-Transform

To enable support for kernel computation via the (improved) Fast-Gauss-Transform
we need the ```figtree``` library. This can be obtained from the original author
Vlad I. Morariu from [here](ttps://github.com/vmorariu/figtree) or you use our
fork (which contains a few small fixes to eliminate compiler warnings) which was
checked-out during the ```git submodule update``` into ```code/figtree```.

Next, compile ```figtree``` to build a *static* library. This is done, since
we will call (in our experiments) ```diagram_distance``` from MATLAB and we
want to avoid having to set the library path. To compile ```figtree``` simply
type

```bash
cd code/figtree
make FIGTREE_LIB_TYPE=static
```

Once this is done, we can configure ```dipha-pss``` to use the library. This is
done by entering the cmake GUI, enabling the ```USE_FGT``` flag and then setting
the correct paths for the ```FIGTREE_INCLUDE``` and ```FIGTREE_LIB``` variable.
In our example, using the command line this would look like (from scratch)

```bash
cd dipha-pss
mkdir build
cd build
cmake .. -DUSE_FGT=ON \
  -DFIGTREE_LIB=../../figtree/lib \
  -DFIGTREE_INCLUDE=../../figtree/include
make
```

This will compile ```diagram_distance``` and enable the ```--use_fgt```
option.

### Compiling DIPHA

You will need to compile DIPHA (i.e., the submodule that we checked out
earlier) in case you want to compute your own persistence diagrams. DIPHA
also uses ```cmake```, so the process should be fairly simple. The
standard process would look like

```bash
cd code/dipha
mkdir build
cd build
cmake ..
make
```
## Examples

In the following, we show some examples which reproduce some of the results
in the *CVPR* and *NIPS* paper. *Note that these examples use MATLAB code
(also contained in the repository).*

### Timing

We start with a simple timing experiment, where we do not use persistence
diagrams computed from data, but simple create random persistence diagrams
for measuring performance of the kernel.

```MATLAB
cd 'code/matlab';
pl_setup
stat = pl_test_timing('/tmp/test', 50:50:200, 20, 0, 1);
disp(stat);
```

This will run ```diagram_distance``` (1) with and (2) without support for FGT
and output timing results (in seconds) for both variants. In particular, we
start with diagrams of 50 random points and increase this to 200 in steps of
50 points. In every run, 20 such diagrams are created, resulting in 20x20 Gram
matrices. The output is written to ```/tmp/test``` which is created in case it
does not exist.

### Averaging PSS feature maps

In this example, we take a large (random) sample of points from a double annulus,
then draw a couple of *small* random samples from that collection, compute persistence
diagrams (from the small samples) and eventually average the corresponding PSS 
feature maps. The full collection of points and three exemplary random samples are 
shown in the figure below:

![Input](https://github.com/rkwitt/persistence-learning/blob/master/common/pss_averaging_input.png "Input")

This is a good example to illustrate how to compute persistence diagrams from distance 
matrices using DIPHA. In particular, for each random sample, the distance matrix simply is 
the *pairwise Euclidean distance* between all the points in each random sample.

The full functionality is implemented in the MATLAB function ```pl_experiment_pss_average.m```.
To produce the results from the *NIPS 2015* paper (see reference above), we additionally provide a
```.mat``` file ```pl_experiment_pss_average_NIPS15.mat``` which you can load
and pass to the script. This sets the configuration (e.g., radii of annuli, seed,
  etc.) we used in the paper. We run the script as follows:

```matlab
cd 'code/matlab';
pl_setup;
out_dir = '/tmp/out';
load ../../data/pl_experiment_pss_average_NIPS15.mat
result = pl_experiment_pss_average(pl_experiment_pss_average_NIPS15, out_dir);
```

This writes all output files to ```/tmp/out``` including (1) the persistence
diagrams, (2) the distance matrices, (3) plots of all the samples, (4) PSS
feature maps and (5) the average PSS feature map. The computed persistence
diagrams are also available as part of the ```result``` structure that is 
returned by the script.
