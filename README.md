# MILAMIN

MILAMIN is a finite element method implementation in native MATLAB that is capable of doing one million degrees of freedom per minute on a modern desktop computer. This includes pre-processing, solving, and post-processing. The MILAMIN strategies and package are applicable to a broad class of problems in Earth science.

MILAMIN is documented in the G-Cubed article: 
["MILAMIN: MATLAB-based finite element method solver for large problems"
by Marcin Dabrowski, Marcin Krotkiewski, and Daniel W. Schmid, doi:10.1029/2007GC001719](https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2007GC001719). 

### Getting Started
Download the Matlab source code as a zip, unpack on your computer, and run one of the examples, e.g. 
```Matlab
thermal2d_test
```
or
```Matlab
mechanical2d_test
```

### Externals
MILAMIN works better with SuiteSparse package and it requires Triangle for the mesh generation.

 * [SuiteSparse by Tim Davis](http://faculty.cse.tamu.edu/davis/suitesparse.html)
 * [Triangle by Jonathan Richard Shewchuk](https://www.cs.cmu.edu/~quake/triangle.html)

 ### Authors

* Marcin Dabrowski
* Marcin Krotkiewski
* Dani Schmid