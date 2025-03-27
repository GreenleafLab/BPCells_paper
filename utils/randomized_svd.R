
# C++ implementation equivalent to qr.Q(qr(mat))
# This should use less intermediate memory, though,
# reducing garbage collection concerns
# Rcpp::sourceCpp(code="
# #include <RcppEigen.h>

# // [[Rcpp::depends(RcppEigen)]]

# // [[Rcpp::export]]
# Eigen::MatrixXd orthogonalize_cpp(Eigen::MatrixXd mat) {
#     Eigen::ColPivHouseholderQR<Eigen::Ref<Eigen::MatrixXd>> qr(mat);
#     return qr.householderQ() * Eigen::MatrixXd::Identity(mat.rows(), qr.rank());
# }
# ")

Rcpp::sourceCpp(code="
#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

// [[Rcpp::export]]
void orthogonalize_in_place_cpp(Eigen::Map<Eigen::MatrixXd> mat) {
    Rprintf(\"Starting QR decomposition\\n\");
    Eigen::ColPivHouseholderQR<Eigen::Ref<Eigen::Map<Eigen::MatrixXd>>> qr(mat);
    Rprintf(\"Starting Q reconstruction\\n\");
    mat = qr.householderQ() * Eigen::MatrixXd::Identity(mat.rows(), qr.rank());
    Rprintf(\"Finished Q reconstruction\\n\");
}
")

gc_with_logging <- function() {
    cat("Memory usage before gc:\n")
    system("free --human --si")
    gc(verbose=FALSE)
    cat("Memory usage after gc:\n")
    system("free --human --si")
}

#' Randomized subspace construction.
#' Halko et al. Algorithm 4.4 (https://arxiv.org/pdf/0909.4061.pdf)
#' @param A Input matrix (m x n)
#' @param l Columns of output matrix (will be m x l)
#' @param q Iteration count
#' @return Matrix Q (m x l) with orthonormal columns approximating
#'   the range of A
halko_4.4 <- function(A, l, q, frequent_gc=TRUE) {
    m <- nrow(A)
    n <- ncol(A)

    # Due to the need to minimize the number of large matrix copies,
    # all variables Y, Q, and their tilde (~) variants have been
    # collapsed into a single variable name. This makes it easier
    # to have gc() remove any unnecessary variables between steps.
    
    M <- matrix(rnorm(n * l), nrow=n, ncol=l) # 1. Initialize omega (n x l)
    if (frequent_gc) {cat(sprintf("At spot 1 %s\n", Sys.time())); gc_with_logging()}
    M <- A %*% M # 2. Calculate Y0 m x l
    if (frequent_gc) {cat(sprintf("At spot 2 %s\n", Sys.time())); gc_with_logging()}
    orthogonalize_in_place_cpp(M)
    # M <- qr(M) # 2. Compute QR factorization of Y0 to get Q0 and R0
    # if (frequent_gc) {cat(sprintf("At spot 2b %s\n", Sys.time())); gc_with_logging()}
    # M <- qr.Q(M) # 2. Get Q0 (m x l)
    for (j in seq_len(q)) {
        if (frequent_gc) {cat(sprintf("At spot 3, j=%d %s\n", j, Sys.time())); gc_with_logging()}
        M <- t(A) %*% M # 4. Compute ~Y_j (n x l)
        if (frequent_gc) {cat(sprintf("At spot 4, j=%d %s\n", j, Sys.time())); gc_with_logging()}
        orthogonalize_in_place_cpp(M)
        # M <- qr(M) # 4. Compute ~Q_j and ~R_j (n x l)
        # if (frequent_gc) {cat(sprintf("At spot 4b, j=%d %s\n", j, Sys.time())); gc_with_logging()}
        # M <- qr.Q(M) # 4. Get ~Q_j (n x l)
        if (frequent_gc) {cat(sprintf("At spot 5, j=%d  %s\n", j, Sys.time())); gc_with_logging()}
        M <- A %*% M # 5. Compute Y_j (m x l)
        if (frequent_gc) {cat(sprintf("At spot 6, j=%d  %s\n", j, Sys.time())); gc_with_logging()}
        orthogonalize_in_place_cpp(M)
        # M <- qr(M) # 5. Compute Q_j R_j (m x l)
        # if (frequent_gc) {cat(sprintf("At spot 6b, j=%d  %s\n", j, Sys.time())); gc_with_logging()}
        # M <- qr.Q(M)
    }
    M
}

#' Randomized SVD
#'
#' Halko et al. Algorithm 5.1 (https://arxiv.org/pdf/0909.4061.pdf)
#'
#' This algorithm approximates an SVD calculation by performing 2*q + 2
#' dense matrix multiplies, 2*q + 1 QR decompositions, and one SVD. The
#' dense matrix multiplies are all with an l x (m or n) matrix. The QR
#' decompositions are on matrices of m x l (1 + q times) and n x l (q times).
#' The SVD is with a m x l matrix. If A is a matrix such that m > n, we internally
#' transpose to minimize the work required.
#' 
#' @param A Input matrix (m x n)
#' @param k Number of singular vectors to estimate
#' @param l Number of singular vectors to use in intermediate calculations.
#'      (Halko et al. recommend l of k + 5, up to 2*k)
#' @param q Number of power iterations to perform
#' @return List of matrices u (m x k), d (k), and v (n x k) such that
#'    u %*% diag(d) %*% t(v) is approximately equal to A
rsvd <- function(A, k, l=k+10, q=2, frequent_gc=TRUE) {
    if (frequent_gc) {cat(sprintf("At spot 7 %s\n", Sys.time())); gc_with_logging()}
    if (nrow(A) > ncol(A)) {
        # Intermediate matrix is smaller when A has more cols than rows,
        # so use recursive call to handle transpose
        res_t <- rsvd(t(A), k, l, q)
        if (frequent_gc) {cat(sprintf("At spot 8 %s\n", Sys.time())); gc_with_logging()}
        return(list(
            d = res_t$d,
            u = res_t$v,
            v = res_t$u,
            nops=res_t$nops
        ))
    }

    Q <- halko_4.4(A, l, q)
    if (frequent_gc) {cat(sprintf("At spot 9 %s\n", Sys.time())); gc_with_logging()}
    B <- t(Q) %*% A
    if (frequent_gc) {cat(sprintf("At spot 10 %s\n", Sys.time())); gc_with_logging()}
    # res <- svd(B, nu=k, nv=k)
    res <- RSpectra::svds(B, k=k)
    res$niter <- NULL
    res$nconv <- NULL
    if (frequent_gc) {rm(B); cat(sprintf("At spot 11 %s\n", Sys.time())); gc_with_logging()}
    res$u <- Q %*% res$u
    if (frequent_gc) {rm(Q); cat(sprintf("At spot 12 %s\n", Sys.time())); gc_with_logging()}
    res$nops <- 2*q + 2
    return(res)
}
