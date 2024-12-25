# BZX
Working repository for porting Fortran BZX code to Python

### src/ : Minimum source directory
Run the code.
```
mkdir check
mkdir geom
gfortran BZX.f90 -o BZX.exe
./BZX.exe
```
You get outputs in check/ and geom/ directories.

### for_check_output/ : Benchmark data of check/ and geom/.
