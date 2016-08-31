# octuplet_analysis #
Analysis code for the MM octuplet + scintillator data

## How to use ##
Take MM data with MMFE8 boards and scintillator data, combined in the format of "src/combined_example.dat". 
Change the input file in "src/octuplet_ana.C" to point to the data file. This will input method will be revised at a later date. 

To use the analysis code:
    cmd_line$> cd octuplet_analysis
    cmd_line$> make 
    cmd_line$> ./octuplet.x

A file named "test.root" should be outputed. 

