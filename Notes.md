# Readme and compiling
- Readme is clear and with nice examples. Project compiles and run with no issues. Shows nice scalability.

- Some errors in the makefiles: the one in for the 1D version is missing the prerequisites. So it does not recompile the code if you change one of its sources or headers. 

- Instead of (or in addition to) `Time gained: xxxx ms` it would be nice to have a % speedup with respect to the number of threads to see how much the code scales. 
And use microseconds instead of milliseconds if you want more precise measures in codes that are indeeed very fast!. In my PC the timing for the 1D test case is 0 ms! If I use microseconds I get a number different from zero (around 12).    ✔

- Are you working on a MAC? I had some problems with a file that has a
name in all uppercase characters, while in the makefile is indicated lowercase... In MAC-OS file name is case insensitive (don't ask me why). Not in any other Unix-type environments.  ✔

# Code

* In the 2D case, why are you including `mpi.h`, if you are not using mpi? By the way, why your compiler does not give you an error?    ✔  
* You often use `const char *` instead of `std::string` C++ strings have been introduced to have safer (and more powerful) strings, compatible with C-style strings. Better use them.
* You often pass `cVector`, `std::vector` and `std::string` by copy instead of by reference. Prefer `const &` to pass by value, particularly for potentially large objects. And, it any case, it does not harm.     ✔
* You should avoid using  `new, malloc, free, delete` in modern C++, rely instead on std containers and if you need raw data buffers employ the method `data()`. Prefer smart pointers to raw pointers when the pointer is "owning" (i,e, it is supposed to handle the resource).
* Nice the use of traits.

