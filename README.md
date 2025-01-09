# Numerical Tasks for Quantum Computing 
## Luca Gawalleck

### Sheet 3
The `python` file contains my implemented code. It can be compiled with

```python3 sheet3/sheet3.py```

It requires the `numpy` library. Beside the required test $\mathrm{CNOT}_{12}\mathrm{H}_1|00\rangle$ 
it shows also different test. The required one is highlighted with hashtags.

The `R` file contains code using the `qsimulatR` package for comparing the results 
with the own implementation. It can be compiled with

```Rscript sheet3/sheet3.r```


### Sheet 4
Compiling with

```Rscript sheet4/sheet4.r```

package `qsimulatR` is required.

```Rplots.pdf``` contains the circuit for task 6 and the histogram for task 8.


### Sheet 6
Compiling with

```Rscript sheet6/sheet6.r```

package `qsimulatR` is required.


### Sheet 7
For task 1 and 2 the `python` file contains the code implementation. It can be compiled with

```python3 sheet7/sheet7.py```

It requires the `numpy` library.

For task 3 the `R` file contains the code using the `qsimulatR` package. It can be compiled with

```Rscript sheet7/sheet7.r```


### Sheet 8
Compiling with

```Rscript sheet8/sheet8.r``` 

The package `qsimulatR` is required. The plots can be found in `Rplots.pdf` or in the
pdf file for the analytical part.

### Sheet 10
Compiling with

```Rscript sheet10/sheet10.r```

The package `qsimulatR` is required. The plots can be found in `Rplots.pdf` or in the
pdf file for the analytical part.


### Sheet 11
Compiling with

```Rscript sheet11/sheet11.r```

The package `qsimulatR` is required. The found r value is printed. 
The code is very slow in this R implementation. The result of the phase estimation algorithm can therefore be found in the file `sheet11/measurement_results.txt`. To safe time and just see the result, one can comment out the phase estimation part and just run the continued fraction algorithm with the imported values. 

!That is works one has to be in the sheet 11 directory or change the path!

Rplots.pdf contains a histogram of the most likely r values according to the phase estimation and the continued fraction algorithm. 
