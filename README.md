# Decomposition of Robust Closed-Loop Dynamic Real-Time Optimization via Generalized Benders Decomposition

In the paper, we proposed the use of Generalized Benders Decomposition (GBD) in the context of Robust Closed-loop Dynamic Real-time Optimization (CLDRTO). 

Unlike traditional DRTO, CLDRTO considers MPC behavior in economic decision-making. To that end, the MPC problem is recast as algebraic equations based on its first-order Karush-Kuhn-Tucker (KKT) optimality conditions, which are added to the CLDRTO constraint set. The robust optimization problem adopts a 2-stage scenario-based representation for stochastic programming. To assess GBD's impact, we apply it in a bioreactor case study, where the problem scenarios are decomposed into subproblems that can be solved in parallel. Additionally, the use of a constrained MPC in this case study requires the addition of the KKT complementarity slackness condition to the CLDRTO constraints. This change can potentially affect GBD's algorithm convergence as convexity cannot be guaranteed; however, experimental results demonstrate that the solution obtained using GBD and the original monolithic problem are consistent in this case. Furthermore, the fact that the decomposition method allows parallel solution of the subproblems significantly reduces computation time when many scenarios are used.

## Description

The paper presents four subcases to evaluate the algorithm's performance:
1. A linearized plant model with a single reactor is used to assess the decomposition algorithm's capacity to reach an optimal solution when applied to a robust CLDRTO with a target tracking objective function and a varying number of equiprobable scenarios.
2. The first subcase is revisited, but instead of using a target tracking objective function at the CLDRTO level, we apply an economic-like function (maximizing product concentration).
3. A subcase with identical reactors in parallel is used to examine the effect of an increased problem size on the algorithm's performance without changing the complexity or convexity of the problem.
4. Lastly, a nonlinear plant model with a target tracking objective function is used to determine the algorithm's effectiveness when the constraint set is nonconvex.

## Installation
The codes were developed to run in a Jupyter Notebook using a Julia kernel. The optimization problems are solved using IPOPT 3.14.12 (for solving subproblems) and CPLEX 22.1.0 (for solving linearized master problem) embedded in the JuMP package v1.15.1. 

Therefore, to use the codes above, one needs to install: Julia, JuMP, IPOPT and CPLEX. Both solvers require manual installation. IPOPT requires installing a Julia Package, whereas CPLEX requiers obtaining a commercial license . 

## License
MIT License

Copyright (c) [2023] [Jose O.A. Matias and Christopher L.E. Swartz]

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions: The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


## Contact
Jose Otavio Assumpcao Matias: assumpcj@mcmaster.ca // Christopher L.E. Swartz: swartzc@mcmaster.ca


