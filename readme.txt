Have a look at the following function:

https://github.com/CFDEMproject/LIGGGHTS-PUBLIC/blob/master/src/fix_check_timestep_gran.cpp#L233

Task 1: Please explain what this function calculates, why is it needed?

This function was extracted into the code below. We have removed some
code related to meshes and adapted some code to make the code compile and output
something.

Task 2: Test the code and fix any obvious bugs you find. Explain what tools you
used to find the bugs.

Task 3: There are two for loops over all particles in the function
"fun". Do you see any way of optimizing these? If so, implement the
optimization such that the code runs as fast as possible (while maintaining
correct results). Please attach the output (timings) you get for your compiled
code, and also indicate what compilation command / flags you used

Bonus task 4: There is a comment saying "// this is not exact...". Can
you identify what this comment refers to? I.e. what is not exact?
