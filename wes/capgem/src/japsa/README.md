### This is a revision to Japsa.

### Original Japsa
Japsa is a Java Package for Sequence Analysis. As the name implies, Japsa was
written primarily in Java. It contains a large number of ready to run programs
as well as a Java API. Please note that from our newest versions, it is required
Java 8 to compile and run from the source code. Prebuilt releases would need
JVM 1.8.0_144 or newer to run properly.

Details of Japsa can be found
in its documentation hosted on [ReadTheDocs](http://japsa.readthedocs.org/en/latest/index.html)

Japsa is released under the accompanying BSD-like license.


### Revisions:
Only code related to CapSim are extracted.


### Install:
To install CapSim, one can go to the folder containing the source code (CSiTE/wes/capgem) and then run make.
A sample command:
make install INSTALL_DIR=./ MXMEM=4000m SERVER=true
If INSTALL_DIR is specified to be a directory other than ‘./’, this directory should be added to environment variable PATH .
