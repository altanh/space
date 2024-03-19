# Parallelism

**3/29/2024.**
Naively parallelizing the cover iteration loop is not sufficient.
For example, if `R(x,y)` has N elements but `R = {(1,1),...,(1,N)}`,
then an iteration over `R(x)` simply can't be parallelized.

> Consider the parallel strategy SuiteSparse:GraphBLAS uses for dense rows in
> masked-SpGEMM. If the rows are too dense (for example if the LHS has only a
> few rows), then the second loop (the loop over columns) is parallelized with
> atomics for accumulation. Maybe we can generalize this. Needs work estimation.

One generic approach is to exploit distributivity.
To compute `R * S`, partition `R`, `S` into disjoint sums `R1 + R2`, `S1 + S2`.
Then `R * S = (R1 + R2) * (S1 + S2) = R1 * S1 + R1 * S2 + R2 * S1 + R2 * S2`.
Each individual join can be executed in parallel.
There is no load balance guarantee.

Is there a more FJ-native parallelization strategy?
