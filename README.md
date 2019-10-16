# Sturm-Liouville-Eigenvalue-Problem
A fast algorithm to solve Sturm-Liouville Eigenvalue problem with certain restrictions on q.

In this code, we solve Sturm Liouville problem that consists of findign all the eigenfunctions y of L which satisfy boundary conditions of an appropriate type. It operates by constructing a non-oscillatory phase function for the differential equation
                                    y′′(t) + (q(t) − λ)y(t) = 0.

We have considered problems with stronger than usual conditions on the coefficient q. These conditions will ensure not only the existence of solutions, but will enable us to calculate the solutions in an extremely efficient fashion. We will assume that q is positive function given on a compact interval [a, b], and smooth on the open interval? (a, b). Moreover, although our approach can be generalized to other boundary conditions, we will consider only Dirichlet and Neumann conditions.

For more details, please refer to my thesis.
(Thesis is still in progress, will add the link shortly.)
