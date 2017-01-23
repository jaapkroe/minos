MINOS
======================

## Obtaining the Code
    git clone https://gitlab.science.ru.nl/jaapkroes/minos.git

## Compiling
    cd minos
    make 

## Using

    ./minos -h
    ./minos file.xyz

## Testing
The test folder contains
- a scaling test as function of system size
- disk read test, evaluating which read method 
  is fastest for reading text files from disk
- fftw tests for 1D and 2D functions
- an example how to use bitpacking

## Coding
The code has following functions
- compute bond distances and angles
- compute tetrahedral angles
- print (number of) neighbors per atom 

The following functions are partially implemented (see other git branches)
- graph search of rings
- k-space normal and height fluctions correlation functions
