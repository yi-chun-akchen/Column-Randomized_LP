Column-Randomized Linear Programs: Performance Guarantees and Applications
-------------------------------

This repository contains all the code used in the numerical experiments in the paper:

> Y.-C. Akchen and V. V. Mišić (2023). [Column-Randomized Linear Programs: Performance Guarantees and Applications](https://ssrn.com/abstract=3656704). Operations Research, to appear.

Citation:
---------

If you use the code and/or data in this repository in your own research, please cite the above paper as follows:

```bibtex
@article{akchen2023column,
	title={Column-Randomized Linear Programs: Performance Guarantees and Applications},
	author={Akchen, Yi-Chun and Mi\v{s}i\'{c}, Velibor V.},
	journal={Operations Research},
	year={2023},
	note={Available at SSRN: \url{https://ssrn.com/abstract=3656704}
}
```

License:
--------

This code is available under the MIT License.

Copyright (C) 2023 Yi-Chun Akchen and Velibor Misic

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


Software Requirements:
----------------------

To run the code, you need to have:
+ Julia version 0.7.*
+ JuMP version 0.18.*
+ Gurobi version 8.0

The code should be compatible with later versions of Gurobi; it is not compatible with newer versions of Julia/JuMP.


Code Structure:
---------------

The code is structured into two directories:

+ `main/`: Contains the files for the numerical experiments reported in the main text. 
  + `cutting_stock_main.jl`: This file includes most functions related to the cutting stock experiments in Section 5 and Section F.
  + `choice_estimation_main.jl`: This file includes most functions related to the non-parametric choice model estimation experiments in Section 6 and Section G.
  + `figure1.jl`,`table1.jl`,`table2.jl` : These files generate instances that produce the results in the correponding figure/table and save the instances as CSV files. One can control number of generated instances for each problem configuration by changing parameter *n_instance*.


+ `e-companion/`: Contains the files for the numerical experiments reported in the electronic companion (i.e., appendices).
    + `tableEC1.jl`,`tableEC2.jl`,`tableEC3.jl`,`tableEC4.jl`,`tableEC5.jl`: These files generate instances that produce the results in the correponding figure/table and save the instances as CSV files. One can control number of generated instances for each problem configuration by changing parameter *n_instance*. Each file may contain additional functions for the corresponding experiments.
    + `figureEC2.jl`: This file generates the numerics for Figure EC2.
    + `secF4_part1.jl`,`secF4_part2.jl`: These two file generates the outcomes in Section F4. In particular, `secF4_part1.jl` focuses on the exact solution approach via column generation. `secF4_part2.jl` discusses the outcomes from the column-randomized LPs and obtain the number of unique columns, the maximum and average incidence of columns, and other quantities discussed in the section.
