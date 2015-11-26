# Persistence Learning

This code implements the kernel(s) for persistence diagrams proposed in the following
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

# Compilation

The core of the code ```diagram_distance``` depends on DIPHA which is included as
a submodule. After you have checked out the repository to your local harddisk
you can checkout the submodule via

```
git update submodule --recursive
```





