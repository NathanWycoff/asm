norm2 <- function(x) x / sqrt(sum((x)^2))

#' Get the distance between subspaces defined as the ranges of A and B
#'
#' @param A A matrix
#' @param B Another matrix with the same number of rows as A.
#' @return A nonnegative scalar given the norm of the difference between the projectors onto the respective subspaces.
#' @export
subspace_dist <- function(A, B) {
    if (class(A) == "numeric") {
        A <- matrix(A, ncol = 1)
    }
    if (class(B) == "numeric") {
        B <- matrix(B, ncol = 1)
    }
    if(!all(dim(A) == dim(B))) stop("Matrices not of same dimension")
  
    m <- nrow(A)
    #Normalize columns
    # Complete A's columns randomly to form a basis for Rn
    A_b <- cbind(A, matrix(rnorm(m*(m-ncol(A))), ncol = m-ncol(A)))
    # An Orthonormal basis for A's complement.
    A_c <- qr.Q(qr(A_b))[,(ncol(A)+1):m]
    # An orthonormal basis for B
    B_b <- qr.Q(qr(B))

    norm(t(A_c) %*% B_b, '2')
}

#' Estimate the Active Subspace of a Function using Finite Difference Gradients
#'
#' Looks between [-1, 1]
#'
#' @param f The function to eval
#' @param r The max dim of the active subspace
#' @param m The dimension of the underlying/embedding space.
#' @param thresh Bound on angle between subspaces.
#' @param scale Scale all gradients to have norm 1?
#' @return A list with sub, the active subspace, sv, the singular values (all m of them), fs, which gives function values, gs, function grads, and X, which gives sampled locations.
#' @export
fd_est_asm <- function(f, r, m, M = NULL, scale = FALSE) {
    if (is.null(M)) {
        M <- ceiling(5 * r * log(m))
        print(paste("Sampling", M, "locations for", M*(M+1), "many function evaluations..."))
    }
    X <- matrix(runif(m*M, -1, 1), ncol = M)#Each column is a point we're going to evaluate the func at.
    grad_res <- lapply(1:M, function(i) fd_grad(X[,i], f))
    grads <- sapply(grad_res, function(i) i$grad)
    fs <- sapply(grad_res, function(i) i$fx)
    if (scale) {
        grads <- apply(grads, 2, function(i) i / sqrt(sum((i)^2)))
    }
    decomp <- svd(grads)

    return(list(sub = decomp$u[,1:r], sv = decomp$d, fs = fs, gs = grads, X = X))
}

#' Change a function's inputs to live in [-1, 1]
#' 
#' Given an m dimensional function whose inputs live in bounded intervals [a1, b1], ..., [am, bm], return a wrapped version of the function whose inputs live in [-1, 1], ..., [-1, 1].
#'
#' @param f The function to wrap, should have a single vector-valued input.
#' @param domain A list of real tuples, indicating the original domain of the function.
#' @return A function wrapping f.
#' @export
domain_to_unit <- function(f, domain) {
    n11_2_ab <- function(x, a, b) {
        a + (b - a) * (x + 1) / 2
    }
    function(x) {
        xt <- sapply(1:length(x), function(i) n11_2_ab(x[i], domain[[i]][1], domain[[i]][2]))
        f(xt)
    }
}

#' Automatically Select an Active Subspace
#' 
#' Select an active subspace via two methods (or both). Method 'log_gap' Examines for an absolute difference in the base b logarithm of the singular values of greater than thresh. Method 'perc_var' picks singular values until alpha % of the total has been accumulated. Method 'both' selects singular values passing both 'log_gap' and 'perc_var' tests.
#' 
#' @param sv A vector of singular values.
#' @param method One of 'both', 'log_gap', or 'perc_var'. See Details.
#' @param b The log base for the 'log_gap' method.
#' @param thresh The threshold for the 'log_gap' method.
#' @param alpha The proportion of the total that must be exceeded for the 'perc_var' method.
#' @return The index of the smallest singular value to be retained, or NULL if none pass the tests.
#' @export
sel_asm <- function(sv, method = 'both', b = exp(1), thresh = 1, alpha = 0.8) {
    if (!method %in% c('both', 'log_gap', 'perc_var')) {
        stop("'method' must be one of 'both', 'log_gap', or 'perc_var'")
    }
    if (class(sv) != 'numeric') {
        stop("'sv' must be a numeric vector of singular values")
    }
    if (alpha < 0 || alpha > 1) {
        stop("alpha should be in [0,1]")
    }
    if (b < 0) {
        stop("Log base should be positive")
    }
    if (thresh < 0) {
        stop("Threshold should be positive")
    }

    aftgap <- which(abs(diff(log(sv, base = b))) <= thresh)
    aftper <- which(cumsum(sv) / sum(sv) >= alpha)

    if (method == 'log_gap') {
        if (length(aftgap) == 0) {
            warning("No suitable spectral gap discovered")
            return(NULL)
        } else {
            return(aftgap[1] - 1)
        }
    } else if (method == 'perc_var') {
        if (length(aftper) == 0) {
            warning("No suitable spectral gap discovered")
            return(NULL)
        } else {
            return(aftper[1])
        }
    } else if (method == 'both') {
        if (length(aftper) == 0 || length(aftgap) == 0) {
            warning("No suitable spectral gap discovered")
            return(NULL)
        } else {
            return(pmin(aftgap[1] - 1, aftper[1]))
        }
    }
}
