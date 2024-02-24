# libqfml-debug

## Installation instruction

```
$ git clone https://github.com/haydenji0731/libqfml-debug.git --recursive
$ mkdir build && cd build
$ cmake ..
$ make
$ cd ../src
$ cmake -DCMAKE_BUILD_TYPE=Release
$ make -j 4
```
## Examples
As of now, the only program that's been implemented in the rooted tree parser. 
```
# inside src directory, replace the input and output filenames accordingly.
$ rtree_parser input.newick > out
```
You should be able to get a variety of information about your tree such as the number of leaves, internal nodes, as well as a visual representation.
