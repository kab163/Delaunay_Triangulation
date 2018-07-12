# Delaunay_Triangulation
DT code, algorithms, research, etc. for KB

This folder contains all the files related to the C code implementation of the Delaunay Triangulation.

kristi.vtk: Shouldn't need to modify this file. This is the output file that is generated. This is a file that you can open up with visit. It just shows the current triangulation and is a way to check your work.

visit_writer.c: Shouldn't need to modify this file. This is necessary for the function which outputs the file above.

DTcode.C: This is the version of the DT code file that should be modified. It can be compiled with g++ (may have to use the "-std=c++11" flag), but make sure the above file is present in the same folder or you will get a compilation error. 

DelTInput.C: This is a version of the latest DT code that now has input file reading capabilities. It reads in the .txt input file and makes a Delaunay Triangulation accordingly. To give this .C file a new input, make sure you change the name of the input file within main to handle it. Right now it's set to "input.txt". That should be the only change necessary.

input.txt: A sample input file. Currently, this input file holds the x-y coordinates of the 48 US capitals (doesn't include Hawaii or Alaska). This input file is usually used for the Traveling Salesman Problem.

Makefile: a simple makefile with the commands to run DTcode.C

To compile:
make (assuming the makefile above is there)

To run:
./a.out      (or whatever executable name you give it)

To commit/push/pull:
Use git commands. If you'd like you can create a new branch to work from to limit any conflicts. I have a copy of this repo that is separate, so don't worry about overwritting code that shouldn't be - there's a backup.
