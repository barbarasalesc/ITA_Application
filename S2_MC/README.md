# About SimS2 code  

## Overview
This repository contains one of the primary analysis scripts developed during my PhD research within the ReD (Recoil Directionality) experiment. The code is designed to simulate and analyze the ionization signal (S2) of a Liquid Argon (LAr) Time Projection Chamber (TPC).

## Key Features

- Monte Carlo Simulation: Generates LAr TPC ionization signals using the Thomas-Imel recombination model coupled with customizable nuclear stopping power models (Ziegler, Lenz-Jensen, or Molière).

- MCMC Fitter: Implements a Markov-Chain Monte Carlo (MCMC) using a Metropolis-Hastings algorithm to fit data using the Thomas-Imel and Lenz-Jensen models.

## Data Usage & Validation
The MCMC algorithm was thoroughly validated on Monte Carlo samples before being applied to fit real experimental data from the ReD 252-Cf campaign. Because the raw experimental data is proprietary to the collaboration, this repository includes a generated MC sample. This sample serves as a fully functional stand-in, allowing users to run the code and demonstrate the MCMC fitting capabilities.


# Running the script
To run this script, ensure you have the following installed:

- ROOT-CERN: version 6.24/04 or newer

- g++: version 8.3.0 or newer

This script is designed to be executed within the ROOT interactive prompt. To run it, open your terminal and execute:
```bash
$ root
root[0] .L SimS2.C
```


## randomhelp()
This is a helper function for the main random() function. It displays the expected arguments, providing basic usage instructions. To run this function simply run on your terminal:
```cpp
randomhelp()
```
And it will print the necessary arguments to use the function random().

## random()
It is the main function of this script. It is a function given by
```cpp
void random(int n, int opt, const char* funcopt, const char* drawopt, int coloropt)
```
with the arguments being: 
- `int n` = number of events to be simulated
- `int opt` = (0,1,2,3) = (no fluctuation, S2 quantum fluctuations, recoil energy fluctuation + S2 quantum fluctuations, recoil energy fluctuation + S2 quantum fluctuations + S2 detection fluctions + selection cuts)
- `const char* funcopt` = (ziegler, moliere, lenz-jensen)
- `const char *drawopt` = (empty, same)
- `int coloropt` = choose any color number
















