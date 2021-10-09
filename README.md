# Differentially-Private-Nonparametric-Hypothesis-Testing

This was my first C++ project post-"hello world." The folder 

> Differentially-Private-Nonparametric-Hypothesis-Testing/DifferentiallyPrivateNonparametricHypthothesisTesting/ 

contains the implementation of the differentially private analogs of the Mann-Whitney, Kruskal-Wallis, and Wilcoxon signed-rank tests from this paper: https://arxiv.org/pdf/1903.09364.pdf .

The functions are compatible with user-provided random number generators. I implement this with C++20 Concepts, so compiling this code requires C++20 support. 

The folder 

> Differentially-Private-Nonparametric-Hypothesis-Testing/DifPrivNonparametricHypthTestingTests/

includes basic tests for these functions. The tests are also examples for how to compute p-values. You can run the tests by running 

> Differentially-Private-Nonparametric-Hypothesis-Testing/Debug/DifPrivNonparametricHypthTestingTests.exe

**For me, writing this code was just an exercise for learning C++ language features. I did not bother to test that my implementations of the algorithms matched the pseudocode from the above paper. If you plan to use this code, you will need to perform these tests yourself.**
