# Abstract

The Steepest Descent (SD), Conjugate Gradient (CG), Minimal Residual (MR) and Restarted GMRES methods were implemented and the convergence behaviors of these methods were studied for a variant of problems. Furthermore, the eigenvalues accuracies were evaluated for diﬀerent conﬁgurations of the Arnoldi method. The results indicate the required properties of each method are upheld (i.e methods requiring symmetric positive deﬁnite (SPD) matrices do not converge for non-SPD), the eigenvalue accuracy for the Arnoldi method increases with increasing iteration, the Ritz values from the Hessenberg matrix are good approximations for the eigenvalues, the non-reorthogonalized Arnoldi method signiﬁcantly degrades orthogonality of resulting matrices but keep the same eigenvalues, and the shifted-inverse Arnoldi iteration converges to a target faster when speciﬁed.

# Repository Contents

The repository contains the implementations of the Arnoldi, Conjugate Gradient (CG), Restarted-GMRES, Minimum Residual (MR) and Steepest Descent (SD) Methods and the scripts to run the computations on the test matrices from the [University of Florida Sparse Matrix Collection] (https://www.cise.ufl.edu/research/sparse/matrices/list_by_id.html)
