================================================================================
Project:    Hidden Markov Model (HMM) Training (Baum-Welch)
Author:     [Your Name Here]
Roll No:    254101020
Course:     [Your Course Name Here]
================================================================================

## 1. Description

This C++ project implements the Baum-Welch algorithm, a standard solution for the HMM "training" or "re-estimation" problem. 

Given an initial HMM model ($\lambda = (A, B, \pi)$) and a sequence of observations ($O$), this program iteratively adjusts the model parameters to best fit the observation data. It calculates the re-estimated state transition matrix ($A$) and observation emission matrix ($B$) that maximize the probability $P(O | \lambda)$.

The program uses the Forward algorithm, Backward algorithm, and scaling to remain numerically stable.


## 2. Files
* `hmm_training.cpp`: The main C++ source code containing the complete Baum-Welch implementation.
* `Initial_Model.txt`: (Input) The initial HMM model parameters (A, B, and $\pi$).
* `HMM_OBSERVATION_SEQUENCE_1.txt`: (Input) The sequence of observations used for training.
* `HMM_AIJ_FINAL.txt`: (Output) The final, trained state transition (A) matrix.
* `HMM_BJK_FINAL.txt`: (Output) The final, trained observation emission (B) matrix.


## 3. How to Compile and Run

### A. Using Visual Studio
1.  Open the `.sln` or `.vcxproj` file.
2.  Ensure `Initial_Model.txt` and `HMM_OBSERVATION_SEQUENCE_1.txt` are in the project's working directory (usually the same folder as the `.cpp` file).
3.  Build and Run the project (Press F5).
4.  The program will print its progress to the console.
5.  Upon completion, the output files `HMM_AIJ_FINAL.txt` and `HMM_BJK_FINAL.txt` will be created in the same directory.

*(Note: The code includes "StdAfx.h" for compatibility with Visual Studio 2010 default projects. If you have disabled precompiled headers, you can safely delete this include line.)*

### B. Using g++ (Linux/MinGW)
1.  (Optional) If compiling with g++, delete the line `#include "StdAfx.h"` from the top of the `hmm_training.cpp` file.
2.  Open a terminal in the project directory.
3.  Compile the code:
    ```bash
    g++ hmm_training.cpp -o hmm_training -std=c++11
    ```
4.  Run the executable:
    ```bash
    ./hmm_training
    ```
5.  The program will print its progress, and the output files will be generated in the current directory.


## 4. Inputs

* **Initial_Model.txt**: Contains the starting $A$, $B$, and $\pi$ matrices.
* **HMM_OBSERVATION_SEQUENCE_1.txt**: A single line of space-separated integers (from 1 to M=32) representing the observed sequence.


## 5. Outputs

* **HMM_AIJ_FINAL.txt**: The re-estimated $A$ (transition) matrix. The format is an $N \times N$ matrix.
* **HMM_BJK_FINAL.txt**: The re-estimated $B$ (emission) matrix. The format is an $N \times M$ matrix.
* **Console Output**: The program will print the log-probability of the observation sequence for each iteration, allowing you to monitor the model's convergence.