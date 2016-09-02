MINOS
======================

## Obtaining the Code
    git clone https://gitlab.science.ru.nl/jaapkroes/minos.git

## Compiling
    cd minos
    make 

## Testing
In the test folder several cases
- disk read test, evaluating read method 
  is fastest for reading text files from disk
- fftw tests for any 1D or 2D function
- boost graph library use (!broken atm!)

## Coding
The code has following functions
- compute bond distances and angles
- compute tetrahedral angles
- print (number of) neighbors per atom 

The following functions are partially implemented (to be done)
- graph search of rings
- k-space normal and height fluctions correlation functions
