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

- **[Compilation](#compilation)**
  - [Compiling support for the Fast-Gauss-Transform](#compiling-support-for-the-fast-gauss-transform)
  - [Compiling DIPHA](#compiling-dipha)
- **[Examples](#examples)**
  - [Timing](#timing)
  - [Averaging PSS feature maps](#averaging-pss-feature-maps)
  - [Simple classification with SVMs](#simple-classification-with-svms)
  - [Using shapes as input data](#using-shapes-as-input-data)

# Compilation

The core of the code is in ```diagram_distance.cpp``` and depends on
DIPHA which is included as a submodule. After you have checked out the
repository via

```bash
git clone https://github.com/rkwitt/persistence-learning.git
```

you can checkout the submodule(s) via

```bash
git submodule init
git submodule update --recursive
```

In case the submodules do not get checked-out properly (or nothing happens),
execute the helper script

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
checked-out during the ```git submodule update``` into ```code/external/figtree```.

Next, compile ```figtree``` to build a *static* library. This is done, since
we will call (in our experiments) ```diagram_distance``` from MATLAB and we
want to avoid having to set the library path. To compile ```figtree``` simply
type

```bash
cd code/external/figtree
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
  -DFIGTREE_LIB=../../external/figtree/lib \
  -DFIGTREE_INCLUDE=../../external/figtree/include
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
cd code/external/dipha
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

### Simple classification with SVMs

In this demonstration, we will create toy data from the annuli as before, however,
this time the objective is to distinguish samples drawn from a single annulus
and samples drawn from a double-annulus based on their persistence diagrams.

First, we create the sample data. We will use ```/tmp/``` as our output
directory.

```matlab
out_dir = '/tmp';
cd code/matlab;
pl_setup;

cnt=1;
for i=1:10
    % create filename
    filename = fullfile(out_dir, sprintf('dmat_%.3d.dipha', cnt));
    % Draw 100 points, center [0, 0], inner radius = 1, outer radius = 2
    points = pl_sample_annulus([0,0], 1, 2, 100, -1)';
    D = squareform(pdist(points));
    save_distance_matrix(D, filename);
    cnt=cnt+1;
end
for i=1:10
    filename = fullfile(out_dir, sprintf('dmat_%.3d.dipha', cnt));
    % Draw 100 points, two centers [0, 1] and [0, -1]
    points = pl_sample_linked_annuli( ...
        100, [0 1], 1, 1.5, [0 -1], 0.5, 1, -1);
    D = squareform(pdist(points));
    save_distance_matrix(D, filename);
    cnt=cnt+1;
end
```

Next, we compute persistence diagrams (by hand), using DIPHA. The inputs are
(as before) the distance matrices we just created from the point samples. We
use a simple bash script (e.g., ```compute.sh```) to do this job. Just set the variable
```DIPHA_BINARY``` to the correct path to the ```dipha``` binary.

```bash
#!/bin/bash
DIPHA_BINARY=<ADD PATH TO DIPHA BINARY HERE>
for i in `seq 1 20`; do
    SRC_FILE=`printf "dmat_%.3d.dipha" ${i}`
    DST_FILE=`printf "dmat_%.3d.pd" ${i}`
    CMD="${DIPHA_BINARY} --upper_dim 2 ${SRC_FILE} ${DST_FILE}"
    ${CMD}
done
```
Move this script to ```/tmp/``` (since this was the output directory in the
MATLAB code) and execute it:

```bash
cd /tmp/
chmod +x compute.sh
./compute.#!/bin/sh
```
We can now use the PSS kernel to compute the Gram matrix, i.e., the matrix
of pairwise kernel evaluations between all created persistence diagrams.
This can be done very easily, since ```diagram_distance``` accepts a ASCII
file as input, where all persistence diagrams are listed (one per line).
In ```/tmp/``` we create such a list via

```bash
cd /tmp/
find . -name 'dmat*.pd' > diagrams.list
```
Finally, we execute ```diagram_distance``` and compute features up
to dimension two (set via ```---dim```). In our example, we compute
the kernel for 1-dimensional features and set the
time parameter of the kernel (set via ```--time```) to 0.1.


```bash
cd code/dipha-pss/build/bin
./diagram_distance --inner_product --time 0.1 --dim 1 /tmp/diagrams.list > /tmp/kernel.txt
```

The kernel matrix, saved as ```/tmp/kernel.txt``` can now be used, e.g., to
train a SVM classifier. We will use libsvm for that purpose, in particular,
the MATLAB interface to libsvm (see the libsvm documentation on how to compile
the MATLAB interface).

```matlab
labels = [ones(10,1);ones(10,1)*2]; % Create labels for training
pos = randsample(1:20,15); % Indices of diagrams used for training
neg = setdiff(1:20,pos);   % Indices of diagrams used for testing
model = svmtrain(labels(pos),[(1:length(pos))' kernel(pos,pos)], '-t 4 -c 1');
[pred,acc,~] = svmpredict(labels(neg), [(1:5)' kernel(neg,pos)], model);
disp(acc);
```
Ideally, we get an accuracy of 100%, simply because the problem is also
very easy. In the demonstrations that follows, we will use more realistic
data. This example just illustrates the *basic* pipeline when we want to
use the kernel in a classification setup.

### Using shapes as input data

In both the *CVPR* and the *NIPS* paper, we experiment with persistence
diagrams obtained from surfaces of 3D shapes. In particular, filtrations
are computed via sublevel sets of a function defined on a simplical
complex (given by the triangulated surface mesh of the 3D shape in that
case).

![ccmesh](https://github.com/rkwitt/persistence-learning/blob/master/common/ccmesh.png "ccmesh")

In the following steps, we demonstrate a full processing pipeline to 
reproduce the results for one dataset of the *NIPS* paper. In particular,
we go from 3D shapes, represented as surface meshes, to persistence diagrams
and then compute a two-sample hypothesis test. 

We provide the full dataset of 3D corpus callosum shapes from the
*NIPS* paper for research purposes. The datasets of the *CVPR* paper
(i.e., SHREC 2014) can be found online - The processing pipeline is
very similar, except that the hypothesis test is replaced by a 
support vector machine for classification (similar to our previous
experiment).

#### Additional 3rd party code

For the full pipeline to work, we will need some additional MATLAB
code which will make our life easier when dealing with meshes. In
particular, we need:

1. [STLRead](http://www.mathworks.com/matlabcentral/fileexchange/22409-stl-file-reader/content/STLRead/html/stldemo.html)
2. [iso2mesh](http://iso2mesh.sourceforge.net/cgi-bin/index.cgi)
3. [(Scale-Invariant) Heat-Kernel Signature](http://cvn.ecp.fr/personnel/iasonas/code/sihks.zip)

The script ```pl_setup.m``` expects these software packages to be available under
```code/external``` .

#### Data

The data that we use are *segmentations* of the [corpus callosum](https://en.wikipedia.org/wiki/Corpus_callosum), 
i.e., a structure in our brain that connects the two hemispheres. These segmentations are binary
masks (in 3D) for which we also have a surface mesh available (i.e., part of the
output of the segmentation process). The data can be found at:

- [Download](https://drive.google.com/file/d/0BxHF82gaPzgSNUlXeFZvRnJ0MEk/view?usp=sharing) raw meshes (STL files)
- [Download](https://drive.google.com/file/d/0BxHF82gaPzgSWmQyTVZPVDFiN1U/view?usp=sharing) preprocessed meshes for MATLAB
- [Download](https://drive.google.com/file/d/0BxHF82gaPzgSUHNaUGNVVVREdzA/view?usp=sharing) meta data (i.e., subject information)

If you are *not* specifically interested in processing the meshes from scratch,
we recommend using the MATLAB data, since all meshes have been checked
already.

#### Processing pipeline

We use the pre-processed meshes in this example. The functions that will be used
within our scripts are:

- ```utilities/pl_mmd.m```
- ```utilities/pl_normalize_kernel.m```
- ```utilities/pl_mesh2dipha.m```
- ```utilities/pl_mesh2hks.m```
- ```experiments/pl_experiment_OASIS_run_dipha.m```
- ```experiments/pl_experiment_OASIS_run_mmd.m```

First, we download the MATLAB data from
the provided link and save the ```.mat``` file ```OASIS_cc.mat```, e.g., at
```/tmp/OASIS_cc.mat```. To compute simplicial complexes, we then use the
MATLAB function ```pl_experiment_OASIS_run_dipha.m``` in the following way:

```matlab
pl_experiment_OASIS_run_dipha('/tmp/OASIS_cc.mat', 'OASIS_cc', 'cc', '/tmp/output')
```

This will process all meshes in the ```.mat``` file and write the simplicial
complexes, as well as the persistence diagrams into ```/tmp/output``` with
all files prefixed by ```cc_```.

Second, we group all subjects according to the desired membership (here:
*demented* vs. *non-demented* at the first visit), compute the kernel and
finally run the kernel two-sample test. The grouping and visit information
is available as meta-data and should be extracted into ```/tmp/output```.
You will also need to configure the two-sample test and provide the options
(in the form of a MATLAB struct) to ```pl_experiment_OASIS_run_mmd.m```.
An exemplary option file (to reproduce the results of the *NIPS* paper is
provided in the ```data``` directory).

```matlab
load ../data/options_pl_experiment_OASIS_run_dipha.mat
[K,pval] = pl_experiment_OASIS_run_mmd('/tmp/OASIS_cc.mat','OASIS_cc', 'cc', ...
  '/tmp/dipha/Group.txt', ...
  '/tmp/dipha/Visit.txt', ...
  options_pl_experiment_OASIS_run_dipha);
```

#### Mesh preprocessing

In case you don't want to use the preprocessed meshes with already computed
Heat-Kernel signatures (e.g., when you want to set the Heat-Kernel signature
times yourself), unpack the raw data, e.g., to ```/tmp/output``` and also
save the meta-data, e.g., at ```/tmp/output/```.

```matlab
subjects = pl_experiment_OASIS_subjects('/tmp/output/Subjects.txt');
```

Next, edit the ```pl_mesh2hks.m``` MATLAB script and set the desired parameters.
Then, execute

```matlab
OASIS_cc = pl_mesh2hks('/tmp/output/raw_OASIS_cc', subjects, 1000);
```

The final parameter 1000 is the scaling of the mesh (to avoid numerical
instabilities). This produces the same result that is already available
in the MATLAB file ```OASIS_cc.mat```.
