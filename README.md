# Genetic-Algorithm

<b> Universiti Tunku Abdul Rahman (UTAR) - Bachelor of Computer Science (Honours) - UCCC2513 MINI PROJECT </b>

<b> Authors: </b> <br>
1\. Brandon Ting En Junn (21ACB01751) </b> <br>
2\. Ling Ji Xiang (21ACB04584) <br>
3\. Loh Chia Heung (23ACB01684) <br>
4\. Yeap Chun Hong (22ACB06352)

<b> Project Title: </b> [Performance Comparison of Genetic Algorithm with Particle Swam Optimisation](/P14G01.pdf)

This is a research-based mini project that compares and evaluates the performance between Evolutionary Computing techniques such as Genetic Algorithm (GA) and Particle Swam Optimisation (PSO) for optimisation problems. Additionally, different combinations of GA operation techniques are experimented to obtain various results.

The following is the 10 benchmark functions used to evaluate their performance:

1. Sphere Function
2. Ackley Function
3. Rastrigin Function
4. Zakharov Function
5. Axis Parallel Hyper-Ellipsoid Function
6. Griewank Function
7. Sum of Different Powers Function
8. Rotated Hyper-Ellipsoid Function
9. Schwefel 2.22 Function
10. Exponential Function

<i> Note: Their referred mathematical implementation can be found [here](/Appendix%20A%20(Updated).docx). </i>

## Setup
### Source Code
The main source code for this project is structured into 2 directories which are:
- `/Assignment`
    - Contains the prototype GA implementation only.
    - Contains the experiment results for GA only.
    - `GA.cpp` is the main source file.
    - `TextToExcel.ipynb` is used to convert the experiment results from text files (.txt) to Excel files (.xlsx).

- `/Project`
    - Contains the full GA and PSO implementations.
    - Contains the experiment results for GA and PSO.
    - `GA.cpp` and `PSO.cpp` are the main source files.
    - `TextToExcel.ipynb` is used to convert the experiment results from text files (.txt) to Excel files (.xlsx).

### Compilation
Initially, this project was compiled manually to conduct the experiments. However, you can use any other Integrated Development Environment (IDE) for C++ to compile this project.

To generate the `.exe` file with manual compilation:
```sh
g++ -o GA.exe GA.cpp
```
```sh
g++ -o PSO.exe PSO.cpp
```

<i> Note: Make sure to install the C++ Compiler (GCC/C++). </i>

## Results
<i>Note: Our experimental results of this project can be found [here](/P14G01.xlsx). </i>

In conclusion, based on our findings:
- Best GA Model is GA49
    - Dynamic Tournament Selection
    - Uniform Crossover
    - Hybrid Comparing Mutation
    - Combined Replacement

- PSO offers better performance than GA
    - PSO is more suitable for continuous optimisation problems
    - GA is more suitable for discrete optimisation problems
