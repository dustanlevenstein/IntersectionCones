# IntersectionCones

I'm going to make the execution of run.py as user-friendly as possible. On the first running, the user will need to enter the filenames where sage and GAP are installed. The GAP package hecke must also be installed. The user will also enter the path where the IntersectionCones directory can be found. It will be assumed that the subdirectory prog with prog.g and sageprog.g are loaded in that directory.

In addition, on every run the user will be asked which hecke cones they want the GAP program to load, as well as which characteristic p cones. The sage program will then divide the cones into blocks and compute the relevant intersections. Note that the characteristic p cones are the cones we are attempting to approximate here, and they are unknown for large symmetric groups.
