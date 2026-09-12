# EDTSurf source

These files are vendored from CueMol revision
`3173d8af62e211dd37b943ee53b3d6a632e6b5d7`,
`src/modules/surface/edtsurf/`.

EDTSurf: Dong Xu and Yang Zhang, *Generating Triangulated Macromolecular
Surfaces by Euclidean Distance Transform*, PLoS ONE 4(12): e8140 (2009),
<https://doi.org/10.1371/journal.pone.0008140>.

The original permission notice is retained in each source file. Local changes
remove CueMol's common header and use matching C allocation functions for
buffers that the upstream implementation resizes with `realloc`. The smoothing
orientation flag is initialized to false: the reference SES path otherwise
reads uninitialized memory. The wrapper includes every supplied element; the
original hydrogen exclusion assumes EDTSurf's element numbering and incorrectly
excludes phosphorus with CueMol's radius table. Voxel ownership is retained.
