This is all the C++ code I wrote for the two population individual based simulator.

The big chunk of code is "simulation.cpp". It includes all the helper functions, plus the main body.
It's easier to import it into R this way, as no problems with headers arise.
"simulation2.cpp" is essentially the same script, with some later functions added at the end, I kept both for consistency with other sub-folders.
UPDATE: fixed a minor bug, "simulation3.cpp" is the final version.

I have also split all the subfunctions in smaller cpp files to make them more digestable.
Those with the "main_body" suffix cannot be compiled independently as they rely on other helper functions.

"R_sample_code.R" shows you how to import and run the simulator in R.
You should download "simulation.cpp" and add the relevant path in the R script.
