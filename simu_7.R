# Silence only glmnet's standard small-CV-sample notice; keep its fallback unchanged.
cv_glmnet <- function(...) {
  withCallingHandlers(glmnet::cv.glmnet(...), warning = function(w) {
    expected <- paste0("Option grouped=FALSE enforced in cv.glmnet, ",
      "since < 3 observations per fold")
    if (identical(conditionMessage(w), expected)) invokeRestart("muffleWarning")
  })
}

# Simulation7: correctly specified mean regressions and misspecified weights.
# Setting(vii), Supplement Table11. d1=100 excludes the intercept in the code.
# Edit the settings below. The target always uses the complete100-replication pool.
N <- c(1000L, 1800L, 2500L, 5000L)
replication_ids <- 1:100
workers <- 30L
output_dir <- "output"

.simu7_method_names <- c("MQR1", "MQR2", "O-TR'", "QR1", "QR2", "TR", "O-TR", "Oracle")
.simu7_method_keys <- c("MQR1", "MQR2", "O-TRp", "QR1", "QR2", "TR", "O-TR", "Oracle")

# Output directories and numerical helpers.
.prepare_run_output <- function (analysis, configuration, replication_ids, workers, output_dir) 
{
    configuration <- gsub("[^a-z0-9_-]", "_", tolower(configuration))
    if (!dir.exists(output_dir)) 
        dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    if (!dir.exists(output_dir)) 
        stop("Cannot create output directory: ", output_dir)
    output_root <- normalizePath(output_dir, mustWork = TRUE)
    parent <- file.path(output_root, analysis, configuration)
    if (!dir.exists(parent)) 
        dir.create(parent, recursive = TRUE, showWarnings = FALSE)
    if (!dir.exists(parent)) 
        stop("Cannot create configuration directory: ", parent)
    stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
    base_name <- paste0("run_", stamp, "_pid", Sys.getpid())
    created <- FALSE
    for (attempt in 0:999) {
        suffix <- if (attempt == 0L) 
            ""
        else sprintf("_%02d", attempt + 1L)
        run_directory <- file.path(parent, paste0(base_name, suffix))
        if (dir.create(run_directory, showWarnings = FALSE)) {
            created <- TRUE
            break
        }
        if (!dir.exists(run_directory)) 
            stop("Cannot create run directory: ", run_directory)
    }
    if (!created) 
        stop("Cannot allocate a unique run directory.")
    run_directory <- normalizePath(run_directory, mustWork = TRUE)
    layout <- list(run_directory = run_directory, log_directory = file.path(run_directory, "logs"), replication_directory = file.path(run_directory, 
        "replications"), results_directory = file.path(run_directory, "results"))
    for (path in layout[-1L]) {
        if (!dir.create(path)) 
            stop("Cannot create output subdirectory: ", path)
    }
    layout$log_file <- file.path(layout$log_directory, "progress.log")
    if (!file.create(layout$log_file)) 
        stop("Cannot create log file: ", layout$log_file)
    header <- paste0("Analysis: ", analysis, " | ", configuration, "\n", "Replications: ", length(replication_ids), 
        " | workers requested: ", workers, "\n", "Run directory: ", layout$run_directory, "\n", "Log file: ", 
        layout$log_file, "\n", "Follow method-by-method progress in the log file.\n")
    cat(header)
    flush.console()
    flush(stdout())
    cat(header, file = layout$log_file, append = TRUE)
    saveRDS(c(list(analysis = analysis, configuration = configuration, replication_ids = as.integer(replication_ids), 
        workers = as.integer(workers), started_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")), layout), 
        file.path(run_directory, "run_manifest.rds"))
    layout
}

ar1_multiply <- function (x, rho = 0.5) 
{
    left <- right <- numeric(length(x))
    left[1] <- x[1]
    right[length(x)] <- x[length(x)]
    if (length(x) > 1L) {
        for (j in 2:length(x)) left[j] <- x[j] + rho * left[j - 1L]
        for (j in (length(x) - 1L):1L) right[j] <- x[j] + rho * right[j + 1L]
    }
    left + right - x
}
atomic_rds <- function (object, path, compress = TRUE) 
{
    tmp <- paste0(path, ".writing_", Sys.getpid())
    if (file.exists(path) || file.exists(tmp)) 
        stop("Refusing to overwrite: ", path)
    saveRDS(object, tmp, compress = compress)
    if (!file.rename(tmp, path)) 
        stop("Cannot publish: ", path)
    invisible(path)
}
compute_shrinkage_cov <- function (residuals, lambda = "auto", ensure_positive_definite = TRUE) 
{
    if (identical(lambda, "auto")) {
        shrink_cov <- cov.shrink(residuals, verbose = FALSE)
    }
    else {
        shrink_cov <- cov.shrink(residuals, lambda = lambda, verbose = FALSE)
    }
    lambda_val <- attr(shrink_cov, "lambda")
    cov_matrix <- as.matrix(shrink_cov)
    class(cov_matrix) <- "matrix"
    cov_matrix <- 0.5 * (cov_matrix + t(cov_matrix))
    if (ensure_positive_definite) {
        eigen_vals <- eigen(cov_matrix, symmetric = TRUE, only.values = TRUE)$values
        if (any(eigen_vals <= 0)) {
            diag(cov_matrix) <- diag(cov_matrix) + 1e-06
        }
    }
    return(list(matrix = cov_matrix, lambda = lambda_val, lambda.estimated = attr(shrink_cov, "lambda.estimated"), 
        diag_values = diag(cov_matrix)))
}
standardize_data <- function (X_train, X_test = NULL) 
{
    mean_X <- colMeans(X_train)
    sd_X <- apply(X_train, 2, sd)
    X_train_scaled <- sweep(sweep(X_train, 2, mean_X, "-"), 2, sd_X, "/")
    if (!is.null(X_test)) {
        X_test_scaled <- sweep(sweep(X_test, 2, mean_X, "-"), 2, sd_X, "/")
        return(list(train = X_train_scaled, test = X_test_scaled))
    }
    else {
        return(list(train = X_train_scaled))
    }
}
el_drm_raw_q <- function (M) 
{
    M <- as.matrix(M)
    if (is.null(colnames(M))) {
        colnames(M) <- paste0("M", seq_len(ncol(M)))
    }
    d <- ncol(M)
    linear_terms <- M
    colnames(linear_terms) <- paste0("L_", colnames(M))
    square_terms <- M^2
    colnames(square_terms) <- paste0("Q_", colnames(M), "_sq")
    if (d >= 2L) {
        pairs <- combn(seq_len(d), 2)
        interaction_terms <- sapply(seq_len(ncol(pairs)), function(k) {
            M[, pairs[1, k]] * M[, pairs[2, k]]
        })
        if (is.null(dim(interaction_terms))) {
            interaction_terms <- matrix(interaction_terms, ncol = 1)
        }
        colnames(interaction_terms) <- apply(pairs, 2, function(idx) {
            paste0("Q_", colnames(M)[idx[1]], "_", colnames(M)[idx[2]])
        })
        cbind(linear_terms, square_terms, interaction_terms)
    }
    else {
        cbind(linear_terms, square_terms)
    }
}
el_drm_q_transform <- function (M, basis_names = NULL) 
{
    raw_q <- el_drm_raw_q(M)
    if (!is.null(basis_names)) {
        if (!all(basis_names %in% colnames(raw_q))) {
            stop("Some selected basis terms are missing from raw q basis.")
        }
        raw_q <- raw_q[, basis_names, drop = FALSE]
    }
    q_center <- colMeans(raw_q)
    q_scale <- apply(raw_q, 2, sd)
    q_scale[!is.finite(q_scale) | q_scale <= 0] <- 1
    list(center = q_center, scale = q_scale, basis_names = colnames(raw_q))
}
el_drm_apply_q <- function (M, q_transform) 
{
    raw_q <- el_drm_raw_q(M)
    raw_q <- raw_q[, q_transform$basis_names, drop = FALSE]
    sweep(sweep(raw_q, 2, q_transform$center, "-"), 2, q_transform$scale, "/")
}
el_drm_design <- function (X) 
{
    X <- as.matrix(X)
    cbind(`(Intercept)` = 1, X)
}
el_drm_unpack_theta <- function (theta, x_dim, q_dim) 
{
    theta0_end <- x_dim * q_dim
    list(theta0 = matrix(theta[seq_len(theta0_end)], nrow = x_dim, ncol = q_dim), theta1 = matrix(theta[theta0_end + 
        seq_len(theta0_end)], nrow = x_dim, ncol = q_dim))
}
el_drm_profile <- function (theta, x_design, q_mat, A, inner_maxit = 60, inner_tol = 1e-05, lambda_maxit = 50, lambda_tol = 1e-07) 
{
    nuisance <- el_drm_solve_nuisance(theta, x_design, q_mat, A, inner_maxit = inner_maxit, inner_tol = inner_tol, 
        lambda_maxit = lambda_maxit, lambda_tol = lambda_tol)
    n_obs <- nrow(x_design)
    diag_eta <- nuisance$eta[cbind(seq_len(n_obs), seq_len(n_obs))]
    loglik <- sum(log(nuisance$p)) + sum(nuisance$alpha) + sum(diag_eta)
    log_weights <- sweep(sweep(nuisance$eta, 1, nuisance$alpha, "+"), 2, log(nuisance$p), "+")
    weights <- exp(pmax(pmin(log_weights, 700), -700))
    expected_q <- weights %*% q_mat
    resid_q <- q_mat - expected_q
    grad0 <- matrix(0, nrow = ncol(x_design), ncol = ncol(q_mat))
    grad1 <- grad0
    idx0 <- which(A == 0)
    idx1 <- which(A == 1)
    if (length(idx0) > 0) {
        grad0 <- crossprod(x_design[idx0, , drop = FALSE], resid_q[idx0, , drop = FALSE])
    }
    if (length(idx1) > 0) {
        grad1 <- crossprod(x_design[idx1, , drop = FALSE], resid_q[idx1, , drop = FALSE])
    }
    list(loglik = loglik, gradient = c(as.vector(grad0), as.vector(grad1)), nuisance = nuisance)
}
el_drm_solve_lambda <- function (base, q_mat, lambda_start, maxit = 50, tol = 1e-07, min_den = 1e-10) 
{
    lambda <- lambda_start
    for (iter in seq_len(maxit)) {
        den <- as.vector(base + q_mat %*% lambda)
        if (any(!is.finite(den)) || any(den <= min_den)) {
            lambda <- lambda * 0.5
            next
        }
        inv_den <- 1/den
        score <- as.vector(crossprod(q_mat, inv_den))
        if (max(abs(score)) < tol) {
            return(lambda)
        }
        jac <- -crossprod(q_mat * as.vector(inv_den^2), q_mat)
        step <- tryCatch(qr.solve(jac, score), error = function(e) MASS::ginv(jac) %*% score)
        step <- as.vector(step)
        step_scale <- 1
        score_norm <- sum(score^2)
        accepted <- FALSE
        for (ls_iter in seq_len(20)) {
            lambda_try <- lambda - step_scale * step
            den_try <- as.vector(base + q_mat %*% lambda_try)
            if (all(is.finite(den_try)) && all(den_try > min_den)) {
                score_try <- as.vector(crossprod(q_mat, 1/den_try))
                if (sum(score_try^2) <= score_norm || step_scale < 1e-04) {
                  lambda <- lambda_try
                  accepted <- TRUE
                  break
                }
            }
            step_scale <- step_scale * 0.5
        }
        if (!accepted) {
            lambda <- lambda - step_scale * step
        }
    }
    lambda
}
el_drm_solve_nuisance <- function (theta, x_design, q_mat, A, inner_maxit = 60, inner_tol = 1e-05, lambda_maxit = 50, lambda_tol = 1e-07) 
{
    n_obs <- nrow(x_design)
    q_dim <- ncol(q_mat)
    x_dim <- ncol(x_design)
    theta_list <- el_drm_unpack_theta(theta, x_dim, q_dim)
    beta0 <- x_design %*% theta_list$theta0
    beta1 <- x_design %*% theta_list$theta1
    beta_obs <- beta0
    beta_obs[A == 1, ] <- beta1[A == 1, ]
    eta <- beta_obs %*% t(q_mat)
    alpha <- rep(0, n_obs)
    lambda <- rep(0, q_dim)
    p_vec <- rep(1/n_obs, n_obs)
    for (iter in seq_len(inner_maxit)) {
        base_log <- matrixStats::colLogSumExps(sweep(eta, 1, alpha, "+"))
        base <- exp(pmin(base_log, 700))
        den <- as.vector(base + q_mat %*% lambda)
        if (any(!is.finite(den)) || any(den <= 0)) {
            stop("EL-DRM denominator became non-positive")
        }
        p_new <- 1/den
        alpha_new <- -matrixStats::rowLogSumExps(sweep(eta, 2, log(p_new), "+"))
        base_log_new <- matrixStats::colLogSumExps(sweep(eta, 1, alpha_new, "+"))
        base_new <- exp(pmin(base_log_new, 700))
        lambda_new <- el_drm_solve_lambda(base_new, q_mat, lambda, maxit = lambda_maxit, tol = lambda_tol)
        den_new <- as.vector(base_new + q_mat %*% lambda_new)
        if (any(!is.finite(den_new)) || any(den_new <= 0)) {
            stop("EL-DRM lambda update produced non-positive denominator")
        }
        p_final <- 1/den_new
        diff <- max(abs(alpha_new - alpha), abs(lambda_new - lambda), abs(p_final - p_vec))
        alpha <- alpha_new
        lambda <- lambda_new
        p_vec <- p_final
        if (diff < inner_tol) {
            break
        }
    }
    list(p = p_vec, alpha = alpha, lambda = lambda, eta = eta, theta0 = theta_list$theta0, theta1 = theta_list$theta1)
}
el_drm_approx_eval_loss <- function (theta, x_eval, q_eval, A_eval, q_ref, p_ref) 
{
    theta_list <- el_drm_unpack_theta(theta, ncol(x_eval), ncol(q_ref))
    log_p <- log(pmax(p_ref, .Machine$double.xmin))
    add_group <- function(idx, theta_k) {
        if (length(idx) == 0) {
            return(0)
        }
        x_g <- x_eval[idx, , drop = FALSE]
        q_g <- q_eval[idx, , drop = FALSE]
        beta_g <- x_g %*% theta_k
        logits <- sweep(beta_g %*% t(q_ref), 2, log_p, "+")
        sum(matrixStats::rowLogSumExps(logits) - rowSums(beta_g * q_g))
    }
    add_group(which(A_eval == 0), theta_list$theta0) + add_group(which(A_eval == 1), theta_list$theta1)
}
el_drm_approx_loss_grad <- function (theta, x_design, q_mat, A, p_mar) 
{
    theta_list <- el_drm_unpack_theta(theta, ncol(x_design), ncol(q_mat))
    log_p <- log(pmax(p_mar, .Machine$double.xmin))
    loss <- 0
    grad0 <- matrix(0, nrow = ncol(x_design), ncol = ncol(q_mat))
    grad1 <- grad0
    add_group <- function(idx, theta_k) {
        if (length(idx) == 0) {
            return(list(loss = 0, grad = matrix(0, nrow = ncol(x_design), ncol = ncol(q_mat))))
        }
        x_g <- x_design[idx, , drop = FALSE]
        q_obs <- q_mat[idx, , drop = FALSE]
        beta_g <- x_g %*% theta_k
        logits <- sweep(beta_g %*% t(q_mat), 2, log_p, "+")
        log_norm <- matrixStats::rowLogSumExps(logits)
        weights <- exp(sweep(logits, 1, log_norm, "-"))
        expected_q <- weights %*% q_mat
        list(loss = sum(log_norm - rowSums(beta_g * q_obs)), grad = crossprod(x_g, expected_q - q_obs))
    }
    group0 <- add_group(which(A == 0), theta_list$theta0)
    group1 <- add_group(which(A == 1), theta_list$theta1)
    loss <- group0$loss + group1$loss
    grad0 <- group0$grad
    grad1 <- group1$grad
    list(loss = loss, gradient = c(as.vector(grad0), as.vector(grad1)))
}
el_drm_fit_marginal_p <- function (q_mat, train_A, marginal_maxit = 80, inner_maxit = 60, inner_tol = 1e-05, lambda_maxit = 50, 
    lambda_tol = 1e-07) 
{
    x_marginal <- matrix(1, nrow = nrow(q_mat), ncol = 1)
    train_A <- as.numeric(train_A)
    q_dim <- ncol(q_mat)
    theta_init <- rep(0, 2 * q_dim)
    last_profile <- NULL
    objective <- function(theta) {
        prof <- tryCatch(el_drm_profile(theta, x_marginal, q_mat, train_A, inner_maxit = inner_maxit, 
            inner_tol = inner_tol, lambda_maxit = lambda_maxit, lambda_tol = lambda_tol), error = function(e) NULL)
        if (is.null(prof) || !is.finite(prof$loglik)) {
            return(1e+100)
        }
        last_profile <<- prof
        -prof$loglik
    }
    gradient <- function(theta) {
        prof <- tryCatch(el_drm_profile(theta, x_marginal, q_mat, train_A, inner_maxit = inner_maxit, 
            inner_tol = inner_tol, lambda_maxit = lambda_maxit, lambda_tol = lambda_tol), error = function(e) NULL)
        if (is.null(prof) || !is.finite(prof$loglik) || any(!is.finite(prof$gradient))) {
            return(rep(0, length(theta)))
        }
        last_profile <<- prof
        -prof$gradient
    }
    opt <- optim(theta_init, objective, gradient, method = "L-BFGS-B", control = list(maxit = marginal_maxit))
    final_profile <- el_drm_profile(opt$par, x_marginal, q_mat, train_A, inner_maxit = inner_maxit, inner_tol = inner_tol, 
        lambda_maxit = lambda_maxit, lambda_tol = lambda_tol)
    p_mar <- final_profile$nuisance$p
    p_mar/sum(p_mar)
}
el_drm_fit_density_ratio <- function (train_X, train_M, train_A, basis_names = NULL, ridge_grid = c(0, 0.01, 0.03, 0.05, 0.1, 0.2, 
    0.3, 0.5, 1), ridge_cv_folds = 3, ridge_cv_seed = 123, eta = NULL, outer_maxit = 80, marginal_maxit = 80, 
    inner_maxit = 60, inner_tol = 1e-05, lambda_maxit = 50, lambda_tol = 1e-07) 
{
    train_X <- as.matrix(train_X)
    train_M <- as.matrix(train_M)
    train_A <- as.numeric(train_A)
    fit_theta <- function(x_fit, q_fit, A_fit, p_fit, ridge_lambda = 0, eta_fit = NULL) {
        n_fit <- nrow(x_fit)
        x_dim_fit <- ncol(x_fit)
        q_dim_fit <- ncol(q_fit)
        if (is.null(eta_fit)) {
            eta_fit <- 1/max(q_dim_fit + ridge_lambda, 1)
        }
        theta_cur <- rep(0, 2 * x_dim_fit * q_dim_fit)
        y_cur <- theta_cur
        t_cur <- 1
        eta_cur <- eta_fit
        objective <- function(theta, prof = NULL) {
            if (is.null(prof)) {
                prof <- el_drm_approx_loss_grad(theta, x_fit, q_fit, A_fit, p_fit)
            }
            if (is.null(prof) || !is.finite(prof$loss)) 
                return(Inf)
            prof$loss/n_fit + 0.5 * ridge_lambda * sum(theta^2)
        }
        prof_cur <- el_drm_approx_loss_grad(theta_cur, x_fit, q_fit, A_fit, p_fit)
        objective_prev <- objective(theta_cur, prof_cur)
        for (iter in seq_len(outer_maxit)) {
            prof_y <- el_drm_approx_loss_grad(y_cur, x_fit, q_fit, A_fit, p_fit)
            objective_y <- objective(y_cur, prof_y)
            if (is.null(prof_y) || !is.finite(objective_y) || any(!is.finite(prof_y$gradient)) || objective_y > 
                objective_prev + 1e-10) {
                y_cur <- theta_cur
                t_cur <- 1
                prof_y <- el_drm_approx_loss_grad(y_cur, x_fit, q_fit, A_fit, p_fit)
                objective_y <- objective(y_cur, prof_y)
                if (is.null(prof_y) || !is.finite(objective_y) || any(!is.finite(prof_y$gradient))) {
                  break
                }
            }
            grad <- prof_y$gradient/n_fit + ridge_lambda * y_cur
            eta_trial <- eta_cur
            accepted <- FALSE
            for (bt in seq_len(20)) {
                theta_trial <- y_cur - eta_trial * grad
                prof_trial <- tryCatch(el_drm_approx_loss_grad(theta_trial, x_fit, q_fit, A_fit, p_fit), 
                  error = function(e) NULL)
                objective_trial <- if (is.null(prof_trial)) {
                  Inf
                }
                else {
                  objective(theta_trial, prof_trial)
                }
                if (is.finite(objective_trial) && objective_trial <= objective_y - 1e-04 * eta_trial * 
                  sum(grad^2)) {
                  accepted <- TRUE
                  break
                }
                eta_trial <- eta_trial/2
            }
            if (!accepted) {
                break
            }
            theta_next <- theta_trial
            objective_next <- objective_trial
            eta_cur <- eta_trial
            if (max(abs(theta_next - theta_cur)) < 1e-05) {
                theta_cur <- theta_next
                objective_prev <- objective_next
                break
            }
            t_next <- (1 + sqrt(1 + 4 * t_cur^2))/2
            y_next <- theta_next + ((t_cur - 1)/t_next) * (theta_next - theta_cur)
            theta_cur <- theta_next
            y_cur <- y_next
            t_cur <- t_next
            objective_prev <- objective_next
        }
        prof_final <- el_drm_approx_loss_grad(theta_cur, x_fit, q_fit, A_fit, p_fit)
        list(theta = theta_cur, objective = objective(theta_cur, prof_final), loss = prof_final$loss/n_fit, 
            ridge_lambda = ridge_lambda, eta = eta_cur)
    }
    fold_id <- make_stratified_fold_id(train_A, K = ridge_cv_folds, seed = ridge_cv_seed)
    K_inner <- max(fold_id)
    cv_loss <- matrix(Inf, nrow = length(ridge_grid), ncol = K_inner)
    for (v in seq_len(K_inner)) {
        idx_val <- which(fold_id == v)
        idx_fit <- which(fold_id != v)
        X_fit <- train_X[idx_fit, , drop = FALSE]
        X_val <- train_X[idx_val, , drop = FALSE]
        M_fit <- train_M[idx_fit, , drop = FALSE]
        M_val <- train_M[idx_val, , drop = FALSE]
        A_fit <- train_A[idx_fit]
        A_val <- train_A[idx_val]
        x_fit <- el_drm_design(X_fit)
        x_val <- el_drm_design(X_val)
        q_transform_v <- el_drm_q_transform(M_fit, basis_names = basis_names)
        q_fit <- el_drm_apply_q(M_fit, q_transform_v)
        q_val <- el_drm_apply_q(M_val, q_transform_v)
        p_fit <- el_drm_fit_marginal_p(q_fit, A_fit, marginal_maxit = marginal_maxit, inner_maxit = inner_maxit, 
            inner_tol = inner_tol, lambda_maxit = lambda_maxit, lambda_tol = lambda_tol)
        for (l in seq_along(ridge_grid)) {
            ridge_l <- ridge_grid[l]
            fit_l <- tryCatch(fit_theta(x_fit = x_fit, q_fit = q_fit, A_fit = A_fit, p_fit = p_fit, ridge_lambda = ridge_l, 
                eta_fit = NULL), error = function(e) NULL)
            if (!is.null(fit_l) && all(is.finite(fit_l$theta))) {
                cv_loss[l, v] <- tryCatch(el_drm_approx_eval_loss(theta = fit_l$theta, x_eval = x_val, 
                  q_eval = q_val, A_eval = A_val, q_ref = q_fit, p_ref = p_fit)/length(idx_val), error = function(e) Inf)
            }
        }
    }
    mean_cv_loss <- rowMeans(cv_loss)
    if (all(!is.finite(mean_cv_loss))) {
        stop("All ridge parameters failed in EL-DRM inner CV.")
    }
    mean_cv_loss[!is.finite(mean_cv_loss)] <- Inf
    best_idx <- which.min(mean_cv_loss)
    best_ridge <- ridge_grid[best_idx]
    x_design <- el_drm_design(train_X)
    q_transform <- el_drm_q_transform(train_M, basis_names = basis_names)
    q_mat <- el_drm_apply_q(train_M, q_transform)
    p_mar <- el_drm_fit_marginal_p(q_mat, train_A, marginal_maxit = marginal_maxit, inner_maxit = inner_maxit, 
        inner_tol = inner_tol, lambda_maxit = lambda_maxit, lambda_tol = lambda_tol)
    fit <- fit_theta(x_fit = x_design, q_fit = q_mat, A_fit = train_A, p_fit = p_mar, ridge_lambda = best_ridge, 
        eta_fit = eta)
    full_train_log_ratio <- el_drm_predict_log_ratio(list(theta = fit$theta, p_mar = p_mar, q_transform = q_transform, 
        q_train = q_mat, x_dim = ncol(x_design), q_dim = ncol(q_mat)), train_X, train_M)
    list(theta = fit$theta, p_mar = p_mar, q_transform = q_transform, q_train = q_mat, x_dim = ncol(x_design), 
        q_dim = ncol(q_mat), basis_names = q_transform$basis_names, eta = fit$eta, objective = fit$objective, 
        ridge_lambda = best_ridge, ridge_cv_loss = mean_cv_loss, ridge_grid = ridge_grid, full_train_log_ratio = full_train_log_ratio)
}
el_drm_predict_log_ratio <- function (fit, query_X, query_M) 
{
    x_query <- el_drm_design(query_X)
    q_query <- el_drm_apply_q(query_M, fit$q_transform)
    theta_list <- el_drm_unpack_theta(fit$theta, fit$x_dim, fit$q_dim)
    beta0_query <- x_query %*% theta_list$theta0
    beta1_query <- x_query %*% theta_list$theta1
    eta0_train <- beta0_query %*% t(fit$q_train)
    eta1_train <- beta1_query %*% t(fit$q_train)
    log_p <- log(pmax(fit$p_mar, .Machine$double.xmin))
    alpha0 <- -matrixStats::rowLogSumExps(sweep(eta0_train, 2, log_p, "+"))
    alpha1 <- -matrixStats::rowLogSumExps(sweep(eta1_train, 2, log_p, "+"))
    alpha0 - alpha1 + rowSums((beta0_query - beta1_query) * q_query)
}
make_stratified_fold_id <- function (A, K = 3, seed = 123) 
{
    A <- as.numeric(A)
    K <- min(K, min(table(A)))
    if (K < 2) 
        stop("Not enough observations for inner CV.")
    set.seed(seed)
    fold_id <- integer(length(A))
    for (a in sort(unique(A))) {
        idx <- sample(which(A == a))
        fold_id[idx] <- rep(seq_len(K), length.out = length(idx))
    }
    fold_id
}
fit_mv_rf_cf <- function (X_train, M_train, K = 3, mtry_grid = NULL, nodesize_grid = c(3, 5, 10, 50), ntree_tune = 50, 
    ntree_final = 50, seed = 123) 
{
    X_train <- as.matrix(X_train)
    M_train <- as.matrix(M_train)
    n <- nrow(X_train)
    p <- ncol(X_train)
    d <- ncol(M_train)
    if (is.null(colnames(X_train))) 
        colnames(X_train) <- paste0("X", seq_len(p))
    if (is.null(colnames(M_train))) 
        colnames(M_train) <- paste0("M", seq_len(d))
    if (is.null(mtry_grid)) {
        mtry_grid <- unique(pmax(1L, pmin(p, round(c(sqrt(p), p/10, p/5, p/3)))))
    }
    nodesize_grid <- unique(pmax(1L, pmin(n - 1L, nodesize_grid)))
    tune_grid <- expand.grid(mtry = mtry_grid, nodesize = nodesize_grid)
    K <- min(K, n)
    if (K < 2) 
        stop("Not enough observations for RF inner CV.")
    has_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    if (has_seed) 
        old_seed <- get(".Random.seed", envir = .GlobalEnv)
    on.exit({
        if (has_seed) {
            assign(".Random.seed", old_seed, envir = .GlobalEnv)
        } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
            rm(".Random.seed", envir = .GlobalEnv)
        }
    }, add = TRUE)
    set.seed(seed)
    fold_id <- sample(rep(seq_len(K), length.out = n))
    formula_multi <- as.formula(paste0("cbind(", paste(colnames(M_train), collapse = ","), ") ~ ", paste(colnames(X_train), 
        collapse = "+")))
    cv_loss <- matrix(NA_real_, nrow = nrow(tune_grid), ncol = K)
    for (v in seq_len(K)) {
        idx_val <- which(fold_id == v)
        idx_fit <- which(fold_id != v)
        X_fit <- X_train[idx_fit, , drop = FALSE]
        X_val <- X_train[idx_val, , drop = FALSE]
        M_fit <- M_train[idx_fit, , drop = FALSE]
        M_val <- M_train[idx_val, , drop = FALSE]
        M_scale <- apply(M_fit, 2, sd)
        M_scale[!is.finite(M_scale) | M_scale <= 0] <- 1
        dat_fit <- data.frame(X_fit, M_fit, check.names = FALSE)
        for (g in seq_len(nrow(tune_grid))) {
            fit_g <- rfsrc(formula_multi, data = dat_fit, ntree = ntree_tune, block.size = d, importance = FALSE, 
                forest = TRUE, splitrule = "mv", nodesize = tune_grid$nodesize[g], mtry = tune_grid$mtry[g], 
                bootstrap = "by.root", seed = seed + 1000 * v)
            pred_g <- predict(fit_g, newdata = data.frame(X_val, check.names = FALSE))
            pred_val <- sapply(seq_len(d), function(j) {
                pred_g$regrOutput[[j]][[1]]
            })
            if (d == 1) 
                pred_val <- matrix(pred_val, ncol = 1)
            err_std <- sweep(pred_val - M_val, 2, M_scale, "/")
            cv_loss[g, v] <- mean(err_std^2)
        }
    }
    mean_cv_loss <- rowMeans(cv_loss)
    best_id <- which.min(mean_cv_loss)
    best_mtry <- tune_grid$mtry[best_id]
    best_nodesize <- tune_grid$nodesize[best_id]
    dat_full <- data.frame(X_train, M_train, check.names = FALSE)
    fit_final <- rfsrc(formula_multi, data = dat_full, ntree = ntree_final, block.size = d, importance = FALSE, 
        forest = TRUE, splitrule = "mv", nodesize = best_nodesize, mtry = best_mtry, bootstrap = "by.root", 
        seed = seed + 99999)
    oob_pred <- get.mv.predicted(fit_final, oob = TRUE)
    colnames(oob_pred) <- colnames(M_train)
    list(fit = fit_final, oob_pred = oob_pred, best_mtry = best_mtry, best_nodesize = best_nodesize, 
        cv_loss = mean_cv_loss, tune_grid = tune_grid)
}
select_different_M_components_ks <- function (X, M, A, max_selected = min(5L, ncol(M)), min_selected = 1L, rel_score_min = 0.1, ntree = 200, 
    nodesize = 5, mtry = NULL, seed = 123) 
{
    X <- as.matrix(X)
    M <- as.matrix(M)
    A <- as.numeric(A)
    n <- nrow(X)
    p <- ncol(X)
    d <- ncol(M)
    if (is.null(colnames(X))) 
        colnames(X) <- paste0("X", seq_len(p))
    if (is.null(colnames(M))) 
        colnames(M) <- paste0("M", seq_len(d))
    idx1 <- which(A == 1)
    idx0 <- which(A == 0)
    if (length(idx1) < 20 || length(idx0) < 20) {
        stop("Too few observations in one treatment arm for mediator screening.")
    }
    max_selected <- min(max_selected, d)
    min_selected <- min(min_selected, max_selected)
    if (is.null(mtry)) {
        mtry <- max(1L, round(sqrt(p)))
    }
    mtry <- min(mtry, p)
    formula_multi <- as.formula(paste0("cbind(", paste(colnames(M), collapse = ","), ") ~ ", paste(colnames(X), 
        collapse = "+")))
    fit_arm <- function(idx, arm_seed) {
        dat <- data.frame(X[idx, , drop = FALSE], M[idx, , drop = FALSE], check.names = FALSE)
        fit <- rfsrc(formula_multi, data = dat, ntree = ntree, block.size = d, importance = FALSE, forest = TRUE, 
            splitrule = "mv", nodesize = min(nodesize, length(idx) - 1L), mtry = mtry, bootstrap = "by.root", 
            seed = arm_seed)
        oob_pred <- get.mv.predicted(fit, oob = TRUE)
        if (d == 1) 
            oob_pred <- matrix(oob_pred, ncol = 1)
        pred <- predict(fit, newdata = data.frame(X, check.names = FALSE))
        pred_all <- sapply(seq_len(d), function(j) {
            pred$regrOutput[[j]][[1]]
        })
        if (d == 1) 
            pred_all <- matrix(pred_all, ncol = 1)
        list(pred_all = pred_all, oob_pred = oob_pred)
    }
    fit1 <- fit_arm(idx1, seed + 1)
    fit0 <- fit_arm(idx0, seed + 2)
    mean_diff2 <- colMeans((fit1$pred_all - fit0$pred_all)^2, na.rm = TRUE)
    resid1 <- M[idx1, , drop = FALSE] - fit1$oob_pred
    resid0 <- M[idx0, , drop = FALSE] - fit0$oob_pred
    var1 <- apply(resid1, 2, var, na.rm = TRUE)
    var0 <- apply(resid0, 2, var, na.rm = TRUE)
    positive_var <- c(var0, var1)
    positive_var <- positive_var[is.finite(positive_var) & positive_var > 0]
    if (length(positive_var) == 0) {
        stop("Unable to estimate residual variances for mediator screening.")
    }
    var_floor <- max(1e-08, 1e-06 * median(positive_var))
    var1[!is.finite(var1) | var1 < var_floor] <- var_floor
    var0[!is.finite(var0) | var0 < var_floor] <- var_floor
    score <- 0.5 * (var1/var0 + var0/var1 - 2 + mean_diff2 * (1/var1 + 1/var0))
    score[!is.finite(score)] <- 0
    ranked <- order(score, decreasing = TRUE)
    ranked <- order(score, decreasing = TRUE)
    finite_scores <- score[is.finite(score)]
    if (length(finite_scores) < 3L) {
        stop("Too few finite mediator screening scores.")
    }
    score_med <- median(finite_scores)
    score_mad <- mad(finite_scores, center = score_med, constant = 1)
    if (!is.finite(score_mad) || score_mad <= 0) {
        score_mad <- IQR(finite_scores)/1.349
    }
    if (!is.finite(score_mad) || score_mad <= 0) {
        score_mad <- max(1e-08, 1e-06 * max(finite_scores))
    }
    noise_threshold <- score_med + 2 * score_mad
    eligible <- ranked[is.finite(score[ranked]) & score[ranked] > noise_threshold]
    if (length(eligible) == 0L) {
        selected <- ranked[1L]
    }
    else {
        selected <- eligible[seq_len(min(max_selected, length(eligible)))]
    }
    if (length(selected) < min_selected) {
        selected <- ranked[seq_len(min(min_selected, length(ranked)))]
    }
    k_select <- length(selected)
    selected
}
select_el_drm_basis <- function (X, M, A, stability_reps = 10L, subsample_frac = 0.8, cv_folds = 5L, alpha = 0.5, stability_cutoff = 0.8, 
    linear_penalty_factor = 0.25, seed = 123) 
{
    X <- as.matrix(X)
    M <- as.matrix(M)
    A <- as.numeric(A)
    full_basis <- el_drm_raw_q(M)
    linear_names <- grep("^L_", colnames(full_basis), value = TRUE)
    quadratic_names <- grep("^Q_", colnames(full_basis), value = TRUE)
    if (length(linear_names) == 0L) {
        stop("No linear mediator basis found.")
    }
    if (length(quadratic_names) == 0L) {
        return(list(selected_basis = linear_names, linear_basis = linear_names, selected_quadratic = character(0), 
            quadratic_frequency = numeric(0), successful_reps = 0L))
    }
    selection_count <- setNames(rep(0L, length(quadratic_names)), quadratic_names)
    idx1 <- which(A == 1)
    idx0 <- which(A == 0)
    if (length(idx1) < 20L || length(idx0) < 20L) {
        stop("Too few observations for EL-DRM basis screening.")
    }
    successful_reps <- 0L
    for (r in seq_len(stability_reps)) {
        set.seed(seed + r)
        sub1 <- sample(idx1, size = floor(subsample_frac * length(idx1)), replace = FALSE)
        sub0 <- sample(idx0, size = floor(subsample_frac * length(idx0)), replace = FALSE)
        idx_sub <- c(sub1, sub0)
        X_sub <- X[idx_sub, , drop = FALSE]
        B_sub <- full_basis[idx_sub, , drop = FALSE]
        A_sub <- A[idx_sub]
        Z_sub <- cbind(X_sub, B_sub)
        penalty_factor <- c(rep(1, ncol(X_sub)), rep(linear_penalty_factor, length(linear_names)), rep(1, 
            length(quadratic_names)))
        fold_id <- make_stratified_fold_id(A_sub, K = cv_folds, seed = seed + 1000L + r)
        fit_cv <- tryCatch(cv_glmnet(x = Z_sub, y = A_sub, family = "binomial", alpha = alpha, foldid = fold_id, 
            penalty.factor = penalty_factor, standardize = TRUE, type.measure = "deviance"), error = function(e) NULL)
        if (is.null(fit_cv)) {
            next
        }
        beta_hat <- tryCatch(coef(fit_cv, s = "lambda.1se"), error = function(e) NULL)
        if (is.null(beta_hat)) {
            next
        }
        beta_vec <- as.vector(beta_hat)
        beta_names <- rownames(beta_hat)
        selected_now <- beta_names[is.finite(beta_vec) & abs(beta_vec) > 1e-08]
        selected_quad_now <- intersect(selected_now, quadratic_names)
        successful_reps <- successful_reps + 1L
        if (length(selected_quad_now) > 0L) {
            selection_count[selected_quad_now] <- selection_count[selected_quad_now] + 1L
        }
    }
    if (successful_reps < max(3L, ceiling(stability_reps/2))) {
        stop("Too few successful EL-DRM basis-screening repetitions: ", successful_reps, "/", stability_reps)
    }
    selection_frequency <- selection_count/successful_reps
    selected_quadratic <- names(selection_frequency)[selection_frequency >= stability_cutoff]
    selected_basis <- c(linear_names, selected_quadratic)
    list(selected_basis = selected_basis, linear_basis = linear_names, selected_quadratic = selected_quadratic, 
        quadratic_frequency = sort(selection_frequency, decreasing = TRUE), successful_reps = successful_reps)
}
generate_num_tree <- function (n) 
{
    num_trees = 50
    return(num_trees)
}
generate_tune_grid <- function (n_train, d, task = 0) 
{
    if (n_train < 1 || d < 1) 
        stop("n_train and d must be positive integers.")
    if (!task %in% c(0, 1)) 
        stop("task must be 0 (regression) or 1 (classification).")
    if (task == 0) {
        mtry0 <- max(1, floor(d/3))
        splitrule_candidates <- c("variance", "extratrees")
        node0 <- max(5, floor(0.03 * n_train))
        node_cap <- 30
    }
    else {
        mtry0 <- max(1, floor(sqrt(d)))
        splitrule_candidates <- c("gini", "extratrees")
        node0 <- max(1, floor(0.01 * n_train))
        node_cap <- 10
    }
    mtry_candidates <- unique(pmax(1, pmin(d, round(c(0.5, 1, 1.5) * mtry0))))
    min_node_size_candidates <- unique(pmax(1, pmin(node_cap, round(c(0.5, 1) * node0))))
    tuneGrid <- expand.grid(mtry = mtry_candidates, min.node.size = min_node_size_candidates, splitrule = splitrule_candidates, 
        stringsAsFactors = FALSE)
    tuneGrid <- unique(tuneGrid)
    rownames(tuneGrid) <- NULL
    tuneGrid
}
lasso_var_selection <- function (X, y, task = 0) 
{
    selected_vars = colnames(X)
    selected_vars <- gsub("^X\\.", "", selected_vars)
    selected_vars <- gsub("^M\\.", "", selected_vars)
    return(list(selected_vars = selected_vars))
}
multi_y_var_selection <- function (X, M) 
{
    M = as.matrix(M)
    if (nrow(X) != nrow(M)) {
        stop("X and M must have the same number of rows")
    }
    y_names <- colnames(M)
    if (is.null(y_names)) {
        y_names <- paste0("M", 1:ncol(M))
    }
    all_selected_vars <- list()
    for (i in seq_along(y_names)) {
        y_current <- M[, i]
        rf <- ranger(y ~ ., data = data.frame(X = X, y = y_current), importance = "impurity", num.trees = 100)
        importance <- importance(rf)
        selected_vars <- names(importance)[importance > mean(importance) + sd(importance)]
        selected_vars <- gsub("^X\\.", "", selected_vars)
        selected_vars <- gsub("^M\\.", "", selected_vars)
        all_selected_vars[[y_names[i]]] <- selected_vars
    }
    union_vars <- unique(unlist(all_selected_vars))
    if (length(union_vars) == 0) {
        union_vars <- colnames(X)
    }
    union_vars <- gsub("^X\\.", "", union_vars)
    union_vars <- gsub("^M\\.", "", union_vars)
    return(list(union_vars = union_vars, selected_by_y = all_selected_vars))
}
normalized_truncate_q <- function (log_q, ref_log_q, lower_prob = 0.01, upper_prob = 0.99) 
{
    log_q <- as.numeric(log_q)
    ref_log_q <- as.numeric(ref_log_q)
    if (length(ref_log_q) < 2L) {
        stop("Too few reference observations for q normalization/truncation.")
    }
    if (any(!is.finite(log_q)) || any(!is.finite(ref_log_q))) {
        stop("Non-finite log density-ratio values.")
    }
    if (lower_prob == 0 && upper_prob == 1) {
        m <- max(ref_log_q)
        log_norm_const <- m + log(mean(exp(ref_log_q - m)))
        q_out <- exp(log_q - log_norm_const)
        if (any(!is.finite(q_out))) {
            stop("Non-finite normalized density-ratio values.")
        }
        attr(q_out, "lower_bound") <- 0
        attr(q_out, "upper_bound") <- Inf
        attr(q_out, "norm_const") <- exp(log_norm_const)
        return(q_out)
    }
    q_raw <- exp(log_q)
    q_ref <- exp(ref_log_q)
    q_bounds <- quantile(q_ref, probs = c(lower_prob, upper_prob), names = FALSE, na.rm = TRUE)
    lower <- q_bounds[1]
    upper <- q_bounds[2]
    q_trunc <- pmax(pmin(q_raw, upper), lower)
    q_ref_trunc <- pmax(pmin(q_ref, upper), lower)
    norm_const <- mean(q_ref_trunc)
    if (!is.finite(norm_const) || norm_const <= 0) {
        stop("Invalid q normalization constant.")
    }
    q_out <- q_trunc/norm_const
    attr(q_out, "lower_bound") <- lower
    attr(q_out, "upper_bound") <- upper
    attr(q_out, "norm_const") <- norm_const
    q_out
}
predict_caret_regression_final <- function (fit, newdata) 
{
    pred <- predict(fit$finalModel, data = newdata)
    if (!is.null(pred$predictions)) {
        as.numeric(pred$predictions)
    }
    else {
        as.numeric(predict(fit, newdata = newdata))
    }
}
one_thread <- function () 
{
    Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1", 
        NUMEXPR_NUM_THREADS = "1", R_RANGER_NUM_THREADS = "1", RF_CORES = "1")
    options(ranger.num.threads = 1L, rf.cores = 1L)
    invisible(NULL)
}
tuning_grids <- function () 
{
    list(rho.range0_default = c(seq(0.01, 2, length.out = 46), 50, 100, 1000, 10000), rho.range1_default = c(seq(0.01, 
        2, length.out = 46), 50, 100, 1000, 10000), rho.range2_default = c(seq(0.01, 11, length.out = 46), 
        50, 100, 1000, 10000), rho.range0_default_CF = c(seq(0.01, 2, length.out = 46), 50, 100, 1000, 
        10000), rho.range1_default_CF = c(seq(0.01, 4, length.out = 46), 50, 100, 1000, 10000), rho.range2_default_CF = c(seq(0.01, 
        10, length.out = 46), 50, 100, 1000, 10000))
}
method_label <- function (x) 
ifelse(x == "O-TRp", "O-TR'", x)
capture_tuning <- function (env) 
{
    values <- mget(c("sel.rho_pi_a", "sel.rho_pi_b", "sel.rho_q", "sel.rho_p", "lambda_mu", "lambda_tau1", 
        "lambda_tau"), envir = env, inherits = FALSE)
    names(values) <- c("rho_pi_a", "rho_pi_b", "rho_q", "rho_p", "lambda_mu", "lambda_tau_n", "lambda_tau")
    values
}
validate_record <- function (z, meta, key, id) 
{
    names_expected <- if (key == "Oracle") 
        "Oracle"
    else c(method_label(key), "Oracle")
    stopifnot(identical(z$replication_id, as.integer(id)), identical(z$case_id, meta$config$case_id), 
        identical(z$N, meta$config$N), identical(z$seed, meta$config$seed), identical(z$theta, meta$theta), 
        identical(z$method_names, names_expected), identical(names(z$pred_theta), names_expected), identical(names(z$pred_sig2), 
            names_expected), all(is.finite(z$pred_theta)), all(is.finite(z$pred_sig2)), all(z$pred_sig2 >= 
            0), identical(rownames(z$nuisance_mse), names_expected), identical(colnames(z$nuisance_mse), 
            c("pi_a", "pi_b", "q", "mu", "tau_n", "tau")), all(is.finite(z$nuisance_mse)), all(z$nuisance_mse >= 
            0), all(z$nuisance_mse["Oracle", ] == 0), identical(z$data_md5, meta$cache_manifest$md5[id]), 
        abs(z$pred_theta["Oracle"] - meta$oracle_theta[id]) < 1e-10, abs(z$pred_sig2["Oracle"] - meta$oracle_variance[id]) < 
            1e-08)
    invisible(TRUE)
}
worker_plan <- function (reps, workers) 
{
    data.frame(replication_id = as.integer(reps), worker_slot = (reps - 1L)%%workers + 1L, wave = (reps - 
        1L)%/%workers + 1L)
}
schedule_fixed_reps <- function (cl, plan, fun, args = list(), on_complete = NULL) 
{
    queues <- lapply(seq_along(cl), function(slot) plan$replication_id[plan$worker_slot == slot])
    position <- integer(length(cl))
    active <- rep(NA_integer_, length(cl))
    completed <- integer()
    submit_next <- function(slot) {
        position[slot] <<- position[slot] + 1L
        if (position[slot] > length(queues[[slot]])) {
            active[slot] <<- NA_integer_
            return(invisible(NULL))
        }
        id <- queues[[slot]][position[slot]]
        active[slot] <<- id
        parallel:::sendCall(cl[[slot]], fun, c(list(id = id, slot = as.integer(slot)), args), tag = id)
        invisible(NULL)
    }
    for (slot in seq_along(cl)) submit_next(slot)
    while (length(completed) < nrow(plan)) {
        result <- parallel:::recvOneResult(cl)
        slot <- result$node
        if (inherits(result$value, "try-error")) 
            stop("Worker ", slot, " failed in rep ", result$tag, ": ", result$value)
        if (is.null(slot) || length(slot) != 1L || is.na(active[slot]) || !identical(as.integer(result$tag), 
            active[slot])) 
            stop("Unexpected worker/replication completion.")
        completed <- c(completed, active[slot])
        submit_next(slot)
        if (!is.null(on_complete)) 
            on_complete(completed)
    }
    invisible(completed)
}

# Per-fold progress; no diagnostic fits or random draws.
fold_audit <- function(method, k, index, mu_rows, taun_rows, tau_rows, env) {
  logger::log_info(paste("rep", env$t, c("MQR1", "MQR2")[method], "fold", k, "complete"))
  NULL
}

# Conditional-density estimation and integration for TR.
run_case26_tr <- function (X, M, A, Y, seed, t, Kmb, n_mc = 50L) 
{
    N <- nrow(X)
    n <- floor(N/Kmb)
    A_fac <- factor(A, levels = c(0, 1), labels = c("Control", "Treatment"))
    pred_pi_a <- pred_pi_b <- pred_q <- pred_mu <- pred_tau1 <- matrix(0, 10, N)
    mtry_pi <- mtry_mu <- minsize_pi <- minsize_mu <- numeric()
    ctrl_pi <- caret::trainControl(method = "cv", number = 5, summaryFunction = caret::twoClassSummary, 
        classProbs = TRUE, allowParallel = FALSE)
    ctrl_mu <- caret::trainControl(method = "cv", number = 5, summaryFunction = caret::defaultSummary, 
        classProbs = FALSE, savePredictions = FALSE, allowParallel = FALSE)
    set.seed(seed)
    for (k in seq_len(Kmb)) {
        index <- seq_len(n) + (k - 1L) * n
        index_nk <- seq_len(N)[-index]
        index_nk_A <- index_nk[A[index_nk] == 1]
        index_nk_A0 <- index_nk[A[index_nk] == 0]
        lasso_pi <- lasso_var_selection(X[index_nk, ], A[index_nk], task = 1)
        selected_vars_pi <- lasso_pi$selected_vars
        pi_rf_TR = caret::train(x = X[index_nk, selected_vars_pi, drop = FALSE], y = A_fac[index_nk], 
            method = "ranger", trControl = ctrl_pi, tuneGrid = generate_tune_grid(length(index_nk), length(selected_vars_pi), 
                1), metric = "ROC", importance = "none", num.trees = generate_num_tree(length(index_nk)))
        prediction_pi_TR = predict(pi_rf_TR$finalModel, data = X[index, selected_vars_pi, drop = FALSE])
        pred_pi_a[6, index] = pmax(pmin(prediction_pi_TR$predictions[, 2], 0.99), 0.01)
        pred_pi_b[6, index] = pred_pi_a[6, index]
        pred_pi_a[7, index] = pred_pi_a[6, index]
        pred_pi_b[7, index] = pred_pi_a[6, index]
        pred_pi_a[9, index] = pred_pi_a[6, index]
        pred_pi_b[9, index] = pred_pi_a[6, index]
        mtry_pi = c(mtry_pi, pi_rf_TR$bestTune$mtry)
        minsize_pi = c(minsize_pi, pi_rf_TR$bestTune$min.node.size)
        pi_rf_TR = butcher(pi_rf_TR)
        rm(pi_rf_TR, prediction_pi_TR)
        gc()
        lasso_mu <- lasso_var_selection(cbind(X[index_nk_A, ], M[index_nk_A, ]), Y[index_nk_A], task = 0)
        selected_vars_mu <- lasso_mu$selected_vars
        mu_rf_TR = caret::train(x = cbind(X[index_nk_A, ], M[index_nk_A, ])[, selected_vars_mu, drop = FALSE], 
            y = Y[index_nk_A], method = "ranger", trControl = ctrl_mu, importance = "none", num.trees = generate_num_tree(length(index_nk_A)), 
            tuneGrid = generate_tune_grid(length(index_nk_A), length(selected_vars_mu), 0))
        pred_mu[6, index] = predict(mu_rf_TR$finalModel, data = cbind(X[index, ], M[index, ])[, selected_vars_mu, 
            drop = FALSE])$predictions
        mtry_mu = c(mtry_mu, mu_rf_TR$bestTune$mtry)
        minsize_mu = c(minsize_mu, mu_rf_TR$bestTune$min.node.size)
        selected_M <- 1:ncol(M)
        ds = length(selected_M)
        x_vars <- multi_y_var_selection(X[index_nk, ], M[index_nk, selected_M, drop = FALSE])$union_vars
        selected_A1 = x_vars
        selected_A0 = x_vars
        rf_mean_A1 <- fit_mv_rf_cf(X_train = X[index_nk_A, selected_A1, drop = FALSE], M_train = M[index_nk_A, 
            selected_M, drop = FALSE], K = 3, ntree_tune = 100, ntree_final = 50, seed = seed + 10000 * 
            t + 100 * k + 1)
        mean_rf_A <- rf_mean_A1$fit
        mean_train_A1 <- rf_mean_A1$oob_pred
        residuals_A1 <- M[index_nk_A, selected_M, drop = FALSE] - mean_train_A1
        rf_mean_A0 <- fit_mv_rf_cf(X_train = X[index_nk_A0, selected_A0, drop = FALSE], M_train = M[index_nk_A0, 
            selected_M, drop = FALSE], K = 3, ntree_tune = 100, ntree_final = 50, seed = seed + 10000 * 
            t + 100 * k + 2)
        mean_rf_A0 <- rf_mean_A0$fit
        mean_train_A0 <- rf_mean_A0$oob_pred
        residuals_A0 <- M[index_nk_A0, selected_M, drop = FALSE] - mean_train_A0
        log_info(paste("TR mean RF fold", k, "A1: mtry =", rf_mean_A1$best_mtry, "nodesize =", rf_mean_A1$best_nodesize, 
            "; A0: mtry =", rf_mean_A0$best_mtry, "nodesize =", rf_mean_A0$best_nodesize))
        residuals_A1 <- sweep(residuals_A1, 2, colMeans(residuals_A1), "-")
        residuals_A0 <- sweep(residuals_A0, 2, colMeans(residuals_A0), "-")
        residuals_pool <- rbind(residuals_A1, residuals_A0)
        cov_A1 <- compute_shrinkage_cov(residuals_A1, lambda = "auto")
        cov_A0 <- compute_shrinkage_cov(residuals_A0, lambda = "auto")
        sig_M_A1 <- cov_A1$matrix
        sig_M_A0 <- cov_A0$matrix
        mean_pred_A0_on_A1 <- predict(mean_rf_A0, newdata = data.frame(X[index_nk_A, selected_A0, drop = FALSE]))
        mean_pred_A0_on_A1 <- sapply(1:length(selected_M), function(i) mean_pred_A0_on_A1$regrOutput[[i]][[1]])
        colnames(mean_pred_A0_on_A1) <- colnames(M)[selected_M]
        log_dens_A1_train <- mvtnorm::dmvnorm(M[index_nk_A, selected_M, drop = FALSE] - mean_train_A1, 
            sigma = sig_M_A1, log = TRUE)
        log_dens_A0_train <- mvtnorm::dmvnorm(M[index_nk_A, selected_M, drop = FALSE] - mean_pred_A0_on_A1, 
            sigma = sig_M_A0, log = TRUE)
        log_q_train <- log_dens_A0_train - log_dens_A1_train
        if (!all(is.finite(log_q_train))) {
            stop("Non-finite training density-ratio prediction.")
        }
        mean_pred <- predict(mean_rf_A, newdata = data.frame(X[index, selected_A1, drop = FALSE]))
        mean_pred_A = sapply(1:ds, function(i) {
            mean_pred$regrOutput[[i]][[1]]
        })
        mean_pred <- predict(mean_rf_A0, newdata = data.frame(X[index, selected_A0, drop = FALSE]))
        mean_pred_A0 = sapply(1:ds, function(i) {
            mean_pred$regrOutput[[i]][[1]]
        })
        colnames(mean_pred_A0) <- colnames(M)[selected_M]
        colnames(mean_pred_A) <- colnames(M)[selected_M]
        log_dens_A1 <- mvtnorm::dmvnorm(M[index, selected_M] - mean_pred_A, sigma = sig_M_A1, log = TRUE)
        log_dens_A0 <- mvtnorm::dmvnorm(M[index, selected_M] - mean_pred_A0, sigma = sig_M_A0, log = TRUE)
        log_q_test_raw <- log_dens_A0 - log_dens_A1
        pred_q[6, index] <- normalized_truncate_q(log_q = log_q_test_raw, ref_log_q = log_q_train)
        if (!all(is.finite(pred_q[6, index]))) {
            stop("Non-finite density ratio after exponentiation.")
        }
        m_vars <- as.numeric(sub("^M", "", grep("^M", selected_vars_mu, value = TRUE)))
        m_names <- colnames(M)[m_vars]
        ds <- length(m_vars)
        if (!all(m_names %in% colnames(mean_pred_A0))) {
            stop("m_vars used by mu must be contained in selected_M used for density estimation.")
        }
        mean_pred_tau1 <- mean_pred_A0[, m_names, drop = FALSE]
        sig_M_tau1 <- sig_M_A0[m_names, m_names, drop = FALSE]
        M_samples <- matrix(0, nrow = length(index) * n_mc, ncol = ds)
        for (i in seq_len(length(index))) {
            M_samples[((i - 1) * n_mc + 1):(i * n_mc), ] <- mvtnorm::rmvnorm(n_mc, mean = mean_pred_tau1[i, 
                ], sigma = sig_M_tau1)
        }
        colnames(M_samples) <- m_names
        pred_data <- do.call(rbind, lapply(seq_len(length(index)), function(i) X[index, , drop = FALSE][rep(i, 
            n_mc), ]))
        pred_data <- cbind(pred_data, M_samples)
        mu_samples <- predict(mu_rf_TR$finalModel, data = pred_data[, selected_vars_mu, drop = FALSE])$predictions
        pred_tau1[6, index] <- tapply(mu_samples, rep(seq_len(length(index)), each = n_mc), mean)
        log_info(paste("rep", t, "TR fold", k, "complete"))
    }
    list(pi = pred_pi_a[6, ], q = pred_q[6, ], mu = pred_mu[6, ], tau_n = pred_tau1[6, ], mtry_pi = mtry_pi, 
        mtry_mu = mtry_mu, minsize_pi = minsize_pi, minsize_mu = minsize_mu)
}

# Parametric and nonparametric nuisance estimation.
fit_parametric_method <- function (data, meta, method, log_file, res_dir) 
{
    stopifnot(identical(data$N, meta$config$N), identical(data$seed, meta$config$seed), identical(data$case_id, 
        meta$config$case_id), data$replication_id %in% seq_len(meta$config$tmax), identical(dim(data$X), 
        c(meta$config$N, meta$config$d1)), identical(dim(data$M), c(meta$config$N, meta$config$d2)))
    list2env(meta$config, envir = environment())
    list2env(meta$parameters, envir = environment())
    list2env(meta$grids, envir = environment())
    theta <- meta$theta
    method_names <- c(method, "Oracle")
    t <- data$replication_id
    n <- floor(N/Kmb)
    num_mu_est <- num_pi_est <- num_q_est <- 10L
    num_p_est <- 6L
    sel.rho_pi_a = c()
    sel.rho_pi_b = c()
    sel.rho_q = c()
    sel.rho_p = c()
    lambda_mu = c()
    lambda_tau = c()
    lambda_tau1 = c()
    mtry_pi = c()
    mtry_p = c()
    mtry_mu = c()
    mtry_tau1 = c()
    mtry_tau = c()
    minsize_pi = c()
    minsize_p = c()
    minsize_mu = c()
    minsize_tau1 = c()
    minsize_tau = c()
    pred_pi_a = matrix(0, nrow = num_pi_est, ncol = N)
    pred_pi_b = matrix(0, nrow = num_pi_est, ncol = N)
    pred_q = matrix(0, nrow = num_q_est, ncol = N)
    pred_p = matrix(0, nrow = num_p_est, ncol = N)
    pred_f1 = matrix(0, nrow = 1, ncol = N)
    pred_f0 = matrix(0, nrow = 1, ncol = N)
    pred_mu = matrix(0, nrow = num_mu_est, ncol = N)
    pred_tau = matrix(0, nrow = num_mu_est, ncol = N)
    pred_tau1 = matrix(0, nrow = num_mu_est, ncol = N)
    pred_theta = matrix(0, nrow = num_mu_est, ncol = N)
    X = data$X
    M = data$M
    colnames(X) = paste0("X", 1:ncol(X))
    colnames(M) = paste0("M", 1:ncol(M))
    A = data$A
    Y = data$Y
    A_fac = factor(A, levels = c(0, 1), labels = c("Control", "Treatment"))
    pi = data$pi
    p = data$p
    q = data$q
    tau = data$tau
    tau1 = data$tau1
    mu = data$mu
    set.seed(seed)
    X_poly = X
    S_poly = cbind(X, M)
    M_poly = M
    pred_pi_a[num_pi_est, ] = pi
    pred_pi_b[num_pi_est, ] = pi
    pred_q[num_q_est, ] = q
    pred_mu[num_mu_est, ] = mu
    pred_tau[num_mu_est, ] = tau
    log_appender(appender_file(log_file))
    log_formatter(formatter_paste)
    log_threshold(INFO)
    fold_audits <- list()
    log_info(paste("Time=", t, "case=", case_id, "started"))
    j <- match(method, c("MQR1", "MQR2"))
    stopifnot(!is.na(j))
    if (j == 1L) 
        for (k in 1:Kmb) {
            log_info(paste("Time=", t, "starting fold", k))
            index = (1:n) + (k - 1) * n
            index_nk = (1:N)[-index]
            index_nk1 = index_nk[1:floor(length(index_nk)/2)]
            index_nk2 = index_nk[-(1:floor(length(index_nk)/2))]
            index_train1 = index_nk[1:floor(length(index_nk)/5)]
            index_nk_2345 = index_nk[-(1:floor(length(index_nk)/5))]
            index_train23 = index_nk_2345[1:floor(length(index_nk_2345)/2)]
            index_train45 = index_nk_2345[-(1:floor(length(index_nk_2345)/2))]
            index_train2 = index_train23[1:floor(length(index_train23)/2)]
            index_train3 = index_train23[-(1:floor(length(index_train23)/2))]
            index_train4 = index_train45[1:floor(length(index_train45)/2)]
            index_train5 = index_train45[-(1:floor(length(index_train45)/2))]
            rm(index_nk_2345, index_train23, index_train45)
            rho.range0 = sqrt(log(d1)/length(index_train1)) * rho.range0_default
            boolFalse = FALSE
            while (boolFalse == FALSE) {
                tryCatch({
                  fit_Log1_1 = glm.regu.cv(x = X_poly[index_train1, ], y = A[index_train1], fold = 5, 
                    loss = "cal", rho.seq = rho.range0, ann.init = FALSE, tune.cut = FALSE)
                  boolFalse = TRUE
                }, error = function(e) {
                }, finally = {
                })
            }
            sel.rho_pi_a = c(sel.rho_pi_a, fit_Log1_1$sel.rho[1]/sqrt(log(d1)/length(index_train1)))
            alpha_a = fit_Log1_1$sel.bet[, 1]
            pred_pi_a[1, index] = plogis(alpha_a[1] + X_poly[index, ] %*% alpha_a[-1])
            rho.range0 = sqrt(log(d1)/length(index_train1)) * rho.range0_default
            boolFalse = FALSE
            while (boolFalse == FALSE) {
                tryCatch({
                  fit_Log1_1 = glm.regu.cv(x = -X_poly[index_train1, ], y = 1 - A[index_train1], fold = 5, 
                    loss = "cal", rho.seq = rho.range0, ann.init = FALSE, tune.cut = FALSE)
                  boolFalse = TRUE
                }, error = function(e) {
                }, finally = {
                })
            }
            sel.rho_pi_b = c(sel.rho_pi_b, fit_Log1_1$sel.rho[1]/sqrt(log(d1)/length(index_train2)))
            alpha_b = fit_Log1_1$sel.bet[, 1]
            alpha_b[1] = -alpha_b[1]
            pred_pi_b[1, index] = plogis(alpha_b[1] + X_poly[index, ] %*% alpha_b[-1])
            S_index2 = S_poly[index_train2, ]
            rho.range1 = sqrt(log(d1 + d2)/length(index_train2)) * rho.range1_default
            boolFalse = FALSE
            while (boolFalse == FALSE) {
                tryCatch({
                  fit_Log2_1 = glm.regu.cv(x = scale(-S_index2), y = A[index_train2], iw = as.numeric(A[index_train2] * 
                    (1 + exp(-alpha_a[1] - X_poly[index_train2, ] %*% alpha_a[-1])) + (1 - A[index_train2]) * 
                    (1 + exp(alpha_b[1] + X_poly[index_train2, ] %*% alpha_b[-1]))), fold = 5, loss = "cal", 
                    rho.seq = rho.range1, ann.init = TRUE, tune.cut = FALSE)
                  boolFalse = TRUE
                }, error = function(e) {
                }, finally = {
                })
            }
            sel.rho_q = c(sel.rho_q, fit_Log2_1$sel.rho[1]/sqrt(log(d1 + d2)/length(index_train2)))
            beta_q = fit_Log2_1$sel.bet[, 1]
            mean_S_index2 = colMeans(-S_index2)
            sd_S_index2 = colSds(-S_index2)
            beta_q_after = rep(0, length(beta_q))
            beta_q_after[-1] = beta_q[-1]/sd_S_index2
            beta_q_after[1] = beta_q[1] - sum((mean_S_index2/sd_S_index2) %*% beta_q[-1])
            beta_q_after[1] = -beta_q_after[1]
            beta_q = beta_q_after
            pred_q[1, index] = exp(beta_q[1] + S_poly[index, ] %*% beta_q[-1])
            fit_lasso_mu = cv_glmnet(x = S_poly[index_train3, ], y = Y[index_train3], weights = A[index_train3] * 
                exp(beta_q[1] + S_poly[index_train3, ] %*% beta_q[-1]) * (1 + exp(-alpha_a[1] - X_poly[index_train3, 
                ] %*% alpha_a[-1])), family = "gaussian", alpha = 1)
            pred_mu[1, index] = predict(fit_lasso_mu, newx = S_poly[index, ], s = "lambda.min")
            lambda_mu = c(lambda_mu, fit_lasso_mu$lambda.min/sqrt(log(d1 + d2)/length(index_train3)))
            index_train4_A0 = index_train4[which(A[index_train4] == 0)]
            pred_mu_new1 = as.numeric(predict(fit_lasso_mu, newx = S_poly[index_train4_A0, , drop = FALSE], 
                s = "lambda.min"))
            tau1_constant_response <- sd(pred_mu_new1) == 0 || length(pred_mu_new1) < 8
            tau1_training_mean <- mean(pred_mu_new1)
            if (tau1_constant_response) {
                pred_tau1[1, index] = rep(tau1_training_mean, length(index))
                lambda_tau1 = c(lambda_tau1, -1)
            }
            else {
                fit_lasso_tau1 = cv_glmnet(x = X_poly[index_train4_A0, ], y = pred_mu_new1, weights = exp(alpha_b[1] + 
                  X_poly[index_train4_A0, ] %*% alpha_b[-1]), family = "gaussian", alpha = 1)
                pred_tau1[1, index] = predict(fit_lasso_tau1, newx = X_poly[index, ], s = "lambda.min")
                lambda_tau1 = c(lambda_tau1, fit_lasso_tau1$lambda.min/sqrt(log(d1)/length(index_train4)))
            }
            pred_mu_new2 = as.numeric(predict(fit_lasso_mu, newx = S_poly[index_train5, , drop = FALSE], 
                s = "lambda.min"))
            if (tau1_constant_response) {
                pred_tau1_new1 = rep(tau1_training_mean, length(index_train5))
            }
            else {
                pred_tau1_new1 = predict(fit_lasso_tau1, newx = X_poly[index_train5, ], s = "lambda.min")
            }
            y_pseudo = pred_tau1_new1 + exp(beta_q[1] + S_poly[index_train5, ] %*% beta_q[-1]) * (Y[index_train5] - 
                pred_mu_new2)
            if (sd(y_pseudo) == 0) {
                pred_tau[1, index] = rep(mean(y_pseudo), length(index))
                lambda_tau = c(lambda_tau, 0)
            }
            else {
                fit_lasso_tau = cv_glmnet(x = X_poly[index_train5, ], y = y_pseudo, weights = A[index_train5] * 
                  exp(-alpha_a[1] - X_poly[index_train5, ] %*% alpha_a[-1]), family = "gaussian", alpha = 1)
                pred_tau[1, index] = predict(fit_lasso_tau, newx = X_poly[index, ], s = "lambda.min")
                lambda_tau = c(lambda_tau, fit_lasso_tau$lambda.min/sqrt(log(d1)/length(index_train4)))
            }
            fold_audits[[length(fold_audits) + 1L]] <- fold_audit(1L, k, index, index_train3, index_train4_A0, 
                index_train5, environment())
        }
    else for (k in 1:Kmb) {
        log_info(paste("Time=", t, "starting fold", k))
        index = (1:n) + (k - 1) * n
        index_nk = (1:N)[-index]
        index_nk_A = index_nk[A[index_nk] == 1]
        index_nk_A0 = index_nk[A[index_nk] == 0]
        rho.range0 = sqrt(log(d1)/length(index_nk)) * rho.range0_default
        boolFalse = FALSE
        while (boolFalse == FALSE) {
            tryCatch({
                fit_Log1_1 = glm.regu.cv(x = X_poly[index_nk, ], y = A[index_nk], fold = 5, loss = "cal", 
                  rho.seq = rho.range0, ann.init = FALSE, tune.cut = FALSE)
                boolFalse = TRUE
            }, error = function(e) {
            }, finally = {
            })
        }
        sel.rho_pi_a = c(sel.rho_pi_a, fit_Log1_1$sel.rho[1]/sqrt(log(d1)/length(index_nk)))
        alpha_a = fit_Log1_1$sel.bet[, 1]
        pred_pi_a[2, index] = plogis(alpha_a[1] + X_poly[index, ] %*% alpha_a[-1])
        rho.range0 = sqrt(log(d1)/length(index_nk)) * rho.range0_default
        boolFalse = FALSE
        while (boolFalse == FALSE) {
            tryCatch({
                fit_Log1_1 = glm.regu.cv(x = -X_poly[index_nk, ], y = 1 - A[index_nk], fold = 5, loss = "cal", 
                  rho.seq = rho.range0, ann.init = FALSE, tune.cut = FALSE)
                boolFalse = TRUE
            }, error = function(e) {
            }, finally = {
            })
        }
        sel.rho_pi_b = c(sel.rho_pi_b, fit_Log1_1$sel.rho[1]/sqrt(log(d1)/length(index_nk)))
        alpha_b = fit_Log1_1$sel.bet[, 1]
        alpha_b[1] = -alpha_b[1]
        pred_pi_b[2, index] = plogis(alpha_b[1] + X_poly[index, ] %*% alpha_b[-1])
        S_index2 = S_poly[index_nk, ]
        rho.range1 = sqrt(log(d1 + d2)/length(index_nk)) * rho.range1_default
        boolFalse = FALSE
        while (boolFalse == FALSE) {
            tryCatch({
                fit_Log2_1 = glm.regu.cv(x = scale(-S_index2), y = A[index_nk], iw = as.numeric(A[index_nk]/plogis(alpha_a[1] + 
                  X_poly[index_nk, ] %*% alpha_a[-1]) + (1 - A[index_nk])/(1 - plogis(alpha_b[1] + X_poly[index_nk, 
                  ] %*% alpha_b[-1]))), fold = 5, loss = "cal", rho.seq = rho.range1, ann.init = TRUE, 
                  tune.cut = FALSE)
                boolFalse = TRUE
            }, error = function(e) {
            }, finally = {
            })
        }
        sel.rho_q = c(sel.rho_q, fit_Log2_1$sel.rho[1]/sqrt(log(d1 + d2)/length(index_nk)))
        beta_q = fit_Log2_1$sel.bet[, 1]
        mean_S_index2 = colMeans(-S_index2)
        sd_S_index2 = colSds(-S_index2)
        beta_q_after = rep(0, length(beta_q))
        beta_q_after[-1] = beta_q[-1]/sd_S_index2
        beta_q_after[1] = beta_q[1] - sum((mean_S_index2/sd_S_index2) %*% beta_q[-1])
        beta_q_after[1] = -beta_q_after[1]
        beta_q[1] = -beta_q[1]
        beta_q = beta_q_after
        pred_q[2, index] = exp(beta_q[1] + S_poly[index, ] %*% beta_q[-1])
        fit_lasso_mu = cv_glmnet(x = S_poly[index_nk_A, ], y = Y[index_nk_A], weights = exp(beta_q[1] + 
            S_poly[index_nk_A, ] %*% beta_q[-1]) * (1 + exp(-alpha_a[1] - X_poly[index_nk_A, ] %*% alpha_a[-1])), 
            family = "gaussian", alpha = 1)
        pred_mu[2, index] = predict(fit_lasso_mu, newx = S_poly[index, ], s = "lambda.min")
        lambda_mu = c(lambda_mu, fit_lasso_mu$lambda.min/sqrt(log(d1 + d2)/length(index_nk)))
        pred_mu_new1 = as.numeric(predict(fit_lasso_mu, newx = S_poly[index_nk_A0, , drop = FALSE], s = "lambda.min"))
        tau1_constant_response <- sd(pred_mu_new1) == 0
        tau1_training_mean <- mean(pred_mu_new1)
        if (tau1_constant_response) {
            pred_tau1[2, index] = rep(tau1_training_mean, length(index))
            lambda_tau1 = c(lambda_tau1, -1)
        }
        else {
            fit_lasso_tau1 = cv_glmnet(x = X_poly[index_nk_A0, ], y = pred_mu_new1, weights = exp(alpha_b[1] + 
                X_poly[index_nk_A0, ] %*% alpha_b[-1]), family = "gaussian", alpha = 1)
            pred_tau1[2, index] = predict(fit_lasso_tau1, newx = X_poly[index, ], s = "lambda.min")
            lambda_tau1 = c(lambda_tau1, fit_lasso_tau1$lambda.min/sqrt(log(d1)/length(index_nk)))
        }
        pred_mu_new1 = as.numeric(predict(fit_lasso_mu, newx = S_poly[index_nk_A, , drop = FALSE], s = "lambda.min"))
        if (tau1_constant_response) {
            pred_tau1_new1 = rep(tau1_training_mean, length(index_nk_A))
        }
        else {
            pred_tau1_new1 = predict(fit_lasso_tau1, newx = X_poly[index_nk_A, ], s = "lambda.min")
        }
        y_pseudo = pred_tau1_new1 + exp(beta_q[1] + S_poly[index_nk_A, ] %*% beta_q[-1]) * (Y[index_nk_A] - 
            pred_mu_new1)
        if (sd(y_pseudo) == 0) {
            pred_tau[2, index] = rep(mean(y_pseudo), length(index))
            lambda_tau = c(lambda_tau, -1)
        }
        else {
            fit_lasso_tau = cv_glmnet(x = X_poly[index_nk_A, ], y = y_pseudo, weights = exp(-alpha_a[1] - 
                X_poly[index_nk_A, ] %*% alpha_a[-1]), family = "gaussian", alpha = 1)
            pred_tau[2, index] = predict(fit_lasso_tau, newx = X_poly[index, ], s = "lambda.min")
            lambda_tau = c(lambda_tau, fit_lasso_tau$lambda.min/sqrt(log(d1)/length(index_nk)))
        }
        fold_audits[[length(fold_audits) + 1L]] <- fold_audit(2L, k, index, index_nk_A, index_nk_A0, 
            index_nk_A, environment())
    }
    oracle_score <- A * q/pi * (Y - mu) + (1 - A)/(1 - pi) * (mu - tau) + tau
    score <- A * pred_q[j, ]/pred_pi_a[j, ] * (Y - pred_mu[j, ]) + (1 - A)/(1 - pred_pi_b[j, ]) * (pred_mu[j, 
        ] - pred_tau1[j, ]) + A/pred_pi_a[j, ] * (pred_tau1[j, ] - pred_tau[j, ]) + pred_tau[j, ]
    scores <- rbind(score, Oracle = as.numeric(oracle_score))
    rownames(scores) <- method_names
    pred_theta <- rowMeans(scores)
    pred_sig2 <- rowMeans((scores - pred_theta)^2)
    nuisance_ids <- j
    nuisance_mse <- t(vapply(nuisance_ids, function(j) c(pi_a = mean((pred_pi_a[j, ] - pi)^2), pi_b = mean((pred_pi_b[j, 
        ] - pi)^2), q = mean((pred_q[j, A == 1] - q[A == 1])^2), mu = mean((pred_mu[j, ] - mu)^2), tau_n = mean((pred_tau1[j, 
        ] - tau1)^2), tau = mean((pred_tau[j, ] - tau)^2)), numeric(6)))
    rownames(nuisance_mse) <- method
    nuisance_mse <- rbind(nuisance_mse, Oracle = rep(0, ncol(nuisance_mse)))
    record <- list(replication_id = as.integer(t), case_id = case_id, N = N, seed = seed, theta = theta, 
        method_names = method_names, pred_theta = pred_theta, pred_sig2 = pred_sig2, nuisance_mse = nuisance_mse, 
        data_md5 = meta$cache_manifest$md5[t])
    validate_record(record, meta, method, t)
    atomic_rds(record, file.path(res_dir, sprintf("rep_%03d.rds", t)))
    logger::log_info(paste("rep", t, method, "complete"))
    invisible(TRUE)
}
fit_comparator_method <- function (t, meta, cache_dir, run_dir, method) 
{
    group <- if (method == "O-TRp") 
        "otrp"
    else if (method == "TR") 
        "tr"
    else "qr"
    one_thread()
    suppressPackageStartupMessages({
        library(glmnet)
        library(RCAL)
        library(matrixStats)
        library(logger)
    })
    if (group %in% c("qr", "tr")) 
        suppressPackageStartupMessages({
            library(caret)
            library(ranger)
            library(randomForestSRC)
            library(butcher)
            library(corpcor)
        })
    list2env(meta$config, envir = environment())
    list2env(meta$parameters, envir = environment())
    list2env(meta$grids, envir = environment())
    theta <- meta$theta
    n <- floor(N/Kmb)
    input_file <- file.path(cache_dir, sprintf("rep_%03d.rds", t))
    if (unname(tools::md5sum(input_file)) != meta$cache_manifest$md5[t]) 
        stop("Data checksum mismatch: ", t)
    data <- readRDS(input_file)
    stopifnot(identical(data$replication_id, as.integer(t)), identical(data$case_id, meta$config$case_id))
    X <- data$X
    M <- data$M
    A <- data$A
    Y <- data$Y
    pi <- data$pi
    q <- data$q
    mu <- data$mu
    tau <- data$tau
    colnames(X) <- paste0("X", seq_len(ncol(X)))
    colnames(M) <- paste0("M", seq_len(ncol(M)))
    X_poly <- X
    S_poly <- cbind(X, M)
    M_poly <- M
    A_fac <- factor(A, levels = c(0, 1), labels = c("Control", "Treatment"))
    pred_pi_a <- pred_pi_b <- pred_q <- pred_mu <- pred_tau <- pred_tau1 <- matrix(0, 10, N)
    pred_p <- matrix(0, 6, N)
    sel.rho_pi_a <- sel.rho_pi_b <- sel.rho_q <- sel.rho_p <- lambda_mu <- lambda_tau1 <- lambda_tau <- numeric()
    mtry_pi <- mtry_p <- mtry_mu <- mtry_tau1 <- mtry_tau <- minsize_pi <- minsize_p <- minsize_mu <- minsize_tau1 <- minsize_tau <- numeric()
    log_file <- meta$log_file
    log_appender(appender_file(log_file))
    log_formatter(formatter_paste)
    log_threshold(INFO)
    cat(format(Sys.time()), sprintf(" | rep %03d | %s started\n", t, group))
    if (group == "otrp") {
        set.seed(seed)
        for (k in 1:Kmb) {
            index = (1:n) + (k - 1) * n
            index_nk = (1:N)[-index]
            index_nk1 = index_nk[1:floor(length(index_nk)/2)]
            index_nk2 = index_nk[-(1:floor(length(index_nk)/2))]
            rho.range0 = sqrt(log(d1)/length(index_nk)) * rho.range0_default
            boolFalse = FALSE
            while (boolFalse == FALSE) {
                tryCatch({
                  fit_Log1_1 = glm.regu.cv(x = X_poly[index_nk, ], y = A[index_nk], fold = 5, loss = "cal", 
                    rho.seq = rho.range0, ann.init = FALSE, tune.cut = FALSE)
                  boolFalse = TRUE
                }, error = function(e) {
                }, finally = {
                })
            }
            sel.rho_pi_a = c(sel.rho_pi_a, fit_Log1_1$sel.rho[1]/sqrt(log(d1)/length(index_nk)))
            sel.rho_pi_b = c(sel.rho_pi_b, fit_Log1_1$sel.rho[1]/sqrt(log(d1)/length(index_nk)))
            alpha_a_miles = fit_Log1_1$sel.bet[, 1]
            pred_pi_a[5, index] = plogis(alpha_a_miles[1] + X_poly[index, ] %*% alpha_a_miles[-1])
            pred_pi_b[5, index] = pred_pi_a[5, index]
            pred_pi_a[8, index] = pred_pi_a[5, index]
            pred_pi_b[8, index] = pred_pi_a[5, index]
            rho.range0 = sqrt(log(d1 + d2)/length(index_nk)) * rho.range0_default
            boolFalse = FALSE
            while (boolFalse == FALSE) {
                tryCatch({
                  fit_Log1_1 = glm.regu.cv(x = S_poly[index_nk, ], y = A[index_nk], fold = 5, loss = "cal", 
                    rho.seq = rho.range0, ann.init = FALSE, tune.cut = FALSE)
                  boolFalse = TRUE
                }, error = function(e) {
                }, finally = {
                })
            }
            sel.rho_p = c(sel.rho_p, fit_Log1_1$sel.rho[1]/sqrt(log(d1 + d2)/length(index_nk)))
            alpha_p = fit_Log1_1$sel.bet[, 1]
            pred_p[3, index] = plogis(alpha_p[1] + S_poly[index, ] %*% alpha_p[-1])
            pred_q[5, index] = (1 - pred_p[3, index])/pred_p[3, index] * pred_pi_a[5, index]/(1 - pred_pi_a[5, 
                index])
            pred_p[5, index] = pred_p[3, index]
            pred_q[8, index] = pred_q[5, index]
            fit_lasso_mu = cv_glmnet(x = S_poly[index_nk1[A[index_nk1] == 1], ], y = Y[index_nk1[A[index_nk1] == 
                1]], family = "gaussian", alpha = 1)
            pred_mu[5, index] = predict(fit_lasso_mu, newx = S_poly[index, ], s = "lambda.min")
            lambda_mu = c(lambda_mu, fit_lasso_mu$lambda.min/sqrt(log(d1 + d2)/length(index_nk1)))
            pred_mu_new1_miles = predict(fit_lasso_mu, newx = S_poly[index_nk2[A[index_nk2] == 0], ], 
                s = "lambda.min")
            if (sd(pred_mu_new1_miles) == 0) {
                pred_tau1[5, index] = rep(mean(pred_mu_new1_miles), length(index))
                lambda_tau1 = c(lambda_tau1, 0)
            }
            else {
                fit_lasso_tau1 = cv_glmnet(x = X_poly[index_nk2[A[index_nk2] == 0], ], y = pred_mu_new1_miles, 
                  family = "gaussian", alpha = 1)
                pred_tau1[5, index] = predict(fit_lasso_tau1, newx = X_poly[index, ], s = "lambda.min")
                lambda_tau1 = c(lambda_tau1, fit_lasso_tau1$lambda.min/sqrt(log(d1)/length(index_nk2)))
            }
        }
        ids <- 5L
        method_names <- c("O-TR'", "Oracle")
        pred_tau[5, ] <- pred_tau1[5, ]
        scores <- rbind(`O-TR'` = A/(1 - pred_pi_a[5, ]) * (1 - pred_p[3, ])/pred_p[3, ] * (Y - pred_mu[5, 
            ]) + (1 - A)/(1 - pred_pi_a[5, ]) * (pred_mu[5, ] - pred_tau1[5, ]) + pred_tau1[5, ])
    }
    else if (group == "tr") {
        tr_fit <- run_case26_tr(X, M, A, Y, seed, t, Kmb, n_mc = 50L)
        pred_pi_a[6, ] <- pred_pi_b[6, ] <- tr_fit$pi
        pred_q[6, ] <- tr_fit$q
        pred_mu[6, ] <- tr_fit$mu
        pred_tau1[6, ] <- pred_tau[6, ] <- tr_fit$tau_n
        mtry_pi <- tr_fit$mtry_pi
        mtry_mu <- tr_fit$mtry_mu
        minsize_pi <- tr_fit$minsize_pi
        minsize_mu <- tr_fit$minsize_mu
        ids <- 6L
        method_names <- c("TR", "Oracle")
        scores <- rbind(TR = A * pred_q[6, ]/pred_pi_a[6, ] * (Y - pred_mu[6, ]) + (1 - A)/(1 - pred_pi_a[6, 
            ]) * (pred_mu[6, ] - pred_tau1[6, ]) + pred_tau1[6, ])
    }
    else {
        ctrl_pi = trainControl(method = "cv", number = 5, summaryFunction = twoClassSummary, classProbs = TRUE)
        ctrl_mu = trainControl(method = "cv", number = 5, summaryFunction = defaultSummary, classProbs = FALSE, 
            savePredictions = FALSE)
        ctrl_pi$allowParallel <- FALSE
        ctrl_mu$allowParallel <- FALSE
        sdr_max_M <- 12L
        sdr_ridge_grid <- c(0, 0.01, 0.03, 0.05, 0.1, 0.2, 0.3, 0.5, 1, 3, 10)
        if (method == "QR1") {
            set.seed(seed)
            for (k in 1:Kmb) {
                index = (1:n) + (k - 1) * n
                index_nk = (1:N)[-index]
                qr1_cut <- floor(seq(0, length(index_nk), length.out = 4))
                index_nk1 <- index_nk[seq_len(qr1_cut[2])]
                index_nk2 <- index_nk[seq.int(qr1_cut[2] + 1L, qr1_cut[3])]
                index_nk3 <- index_nk[seq.int(qr1_cut[3] + 1L, qr1_cut[4])]
                index_nk12 <- c(index_nk1, index_nk2)
                index_nk1_A <- index_nk1[A[index_nk1] == 1]
                index_nk1_A0 <- index_nk1[A[index_nk1] == 0]
                index_nk2_A0 <- index_nk2[A[index_nk2] == 0]
                index_nk3_A <- index_nk3[A[index_nk3] == 1]
                lasso_pi <- lasso_var_selection(X[index_nk, ], A[index_nk], task = 1)
                selected_vars_pi <- lasso_pi$selected_vars
                pi_rf = caret::train(x = X[index_nk, selected_vars_pi, drop = FALSE], y = A_fac[index_nk], 
                  method = "ranger", trControl = ctrl_pi, tuneGrid = generate_tune_grid(length(index_nk), 
                    length(selected_vars_pi), 1), metric = "ROC", importance = "none", num.trees = generate_num_tree(length(index_nk)))
                prediction_pi = predict(pi_rf$finalModel, data = X[index, selected_vars_pi, drop = FALSE])
                pred_pi_a[3, index] = pmax(pmin(prediction_pi$predictions[, 2], 0.99), 0.01)
                pred_pi_b[3, index] = pred_pi_a[3, index]
                mtry_pi = c(mtry_pi, pi_rf$bestTune$mtry)
                minsize_pi = c(minsize_pi, pi_rf$bestTune$min.node.size)
                selected_M <- select_different_M_components_ks(X[index_nk12, ], M[index_nk12, ], A[index_nk12], 
                  max_selected = sdr_max_M, min_selected = 1)
                x_vars <- multi_y_var_selection(X[index_nk12, ], M[index_nk12, selected_M, drop = FALSE])$union_vars
                basis_sel <- select_el_drm_basis(X = X[index_nk12, x_vars, drop = FALSE], M = M[index_nk12, 
                  selected_M, drop = FALSE], A = A[index_nk12], stability_reps = 10, subsample_frac = 0.8, 
                  cv_folds = 5, alpha = 0.5, stability_cutoff = 0.8, seed = seed + 20000 * t + 200 * 
                    k)
                selected_basis <- basis_sel$selected_basis
                log_info(paste("QR1 basis fold", k, "; selected_M =", paste(selected_M, collapse = ","), 
                  "; selected_basis =", paste(selected_basis, collapse = ","), "; successful_reps =", 
                  basis_sel$successful_reps, "; quadratic_freq =", paste(names(basis_sel$quadratic_frequency), 
                    round(basis_sel$quadratic_frequency, 2), sep = ":", collapse = ",")))
                f_el_drm_sdr1 <- el_drm_fit_density_ratio(X[index_nk12, x_vars, drop = FALSE], M[index_nk12, 
                  selected_M, drop = FALSE], A[index_nk12], basis_names = selected_basis, ridge_grid = sdr_ridge_grid, 
                  ridge_cv_folds = 5, ridge_cv_seed = seed + 10000 * t + 100 * k + 2)
                log_ratio_index <- el_drm_predict_log_ratio(f_el_drm_sdr1, X[index, x_vars, drop = FALSE], 
                  M[index, selected_M, drop = FALSE])
                ref_log_q_sdr1 <- f_el_drm_sdr1$full_train_log_ratio[A[index_nk12] == 1]
                pred_q[3, index] <- normalized_truncate_q(log_q = log_ratio_index, ref_log_q = ref_log_q_sdr1)
                log_info(paste("QR1 fold", k, "ridge_lambda =", f_el_drm_sdr1$ridge_lambda, "x_dim =", 
                  f_el_drm_sdr1$x_dim, "q_dim =", f_el_drm_sdr1$q_dim, "theta_dim =", length(f_el_drm_sdr1$theta)))
                lasso_mu <- lasso_var_selection(cbind(X[index_nk1_A, ], M[index_nk1_A, ]), Y[index_nk1_A], 
                  task = 0)
                selected_vars_mu <- lasso_mu$selected_vars
                mu_rf = caret::train(x = cbind(X[index_nk1_A, ], M[index_nk1_A, ])[, selected_vars_mu, 
                  drop = FALSE], y = Y[index_nk1_A], method = "ranger", tuneGrid = generate_tune_grid(length(index_nk1_A), 
                  length(selected_vars_mu), 0), trControl = ctrl_mu, importance = "none", num.trees = generate_num_tree(length(index_nk1_A)))
                pred_mu[3, index] = predict(mu_rf$finalModel, data = cbind(X[index, ], M[index, ])[, 
                  selected_vars_mu, drop = FALSE])$predictions
                mtry_mu = c(mtry_mu, mu_rf$bestTune$mtry)
                minsize_mu = c(minsize_mu, mu_rf$bestTune$min.node.size)
                prediction_mu_1 = predict(mu_rf$finalModel, data = cbind(X[index_nk2_A0, ], M[index_nk2_A0, 
                  ])[, selected_vars_mu, drop = FALSE])
                lasso_tau1 <- lasso_var_selection(X[index_nk2_A0, ], prediction_mu_1$predictions, task = 0)
                selected_vars_tau1 <- lasso_tau1$selected_vars
                tau1_rf = caret::train(x = X[index_nk2_A0, selected_vars_tau1, drop = FALSE], y = prediction_mu_1$predictions, 
                  method = "ranger", tuneGrid = generate_tune_grid(length(index_nk2_A0), length(selected_vars_tau1), 
                    0), trControl = ctrl_mu, importance = "none", num.trees = generate_num_tree(length(index_nk2_A0)))
                mtry_tau1 = c(mtry_tau1, tau1_rf$bestTune$mtry)
                minsize_tau1 = c(minsize_tau1, tau1_rf$bestTune$min.node.size)
                pred_tau1[3, index] = predict(tau1_rf$finalModel, data = X[index, selected_vars_tau1, 
                  drop = FALSE])$predictions
                log_ratio_index_nk3_A <- el_drm_predict_log_ratio(f_el_drm_sdr1, X[index_nk3_A, x_vars, 
                  drop = FALSE], M[index_nk3_A, selected_M, drop = FALSE])
                prediction_q_1 <- normalized_truncate_q(log_q = log_ratio_index_nk3_A, ref_log_q = ref_log_q_sdr1)
                prediction_mu_2 = predict(mu_rf$finalModel, data = cbind(X[index_nk3_A, ], M[index_nk3_A, 
                  ])[, selected_vars_mu, drop = FALSE])
                prediction_tau1_1 = predict(tau1_rf$finalModel, data = X[index_nk3_A, selected_vars_tau1, 
                  drop = FALSE])
                pred_pseudoY = prediction_q_1 * (Y[index_nk3_A] - prediction_mu_2$predictions) + prediction_tau1_1$predictions
                lasso_tau <- lasso_var_selection(X[index_nk3_A, ], pred_pseudoY, task = 0)
                selected_vars_tau <- lasso_tau$selected_vars
                tau_rf = caret::train(x = X[index_nk3_A, selected_vars_tau, drop = FALSE], y = pred_pseudoY, 
                  method = "ranger", tuneGrid = generate_tune_grid(length(index_nk3_A), length(selected_vars_tau), 
                    0), trControl = ctrl_mu, importance = "none", num.trees = generate_num_tree(length(index_nk3_A)))
                mtry_tau = c(mtry_tau, tau_rf$bestTune$mtry)
                minsize_tau = c(minsize_tau, tau_rf$bestTune$min.node.size)
                pred_tau[3, index] = predict(tau_rf$finalModel, data = X[index, selected_vars_tau, drop = FALSE])$predictions
                pi_rf = butcher(pi_rf)
                rm(pi_rf, prediction_pi)
                mu_rf = butcher(mu_rf)
                tau1_rf = butcher(tau1_rf)
                tau_rf = butcher(tau_rf)
                rm(mu_rf, tau1_rf, tau_rf, prediction_mu_1, prediction_mu_2, prediction_tau1_1)
                rm(prediction_q_1)
                gc()
            }
            log_info(paste("rep", t, "QR1 complete"))
            cat(format(Sys.time()), sprintf(" | rep %03d | QR1 complete\n", t))
        }
        if (method == "QR2") {
            set.seed(seed)
            for (k in 1:Kmb) {
                index = (1:n) + (k - 1) * n
                index_nk = (1:N)[-index]
                index_nk_A = index_nk[A[index_nk] == 1]
                index_nk_A0 = index_nk[A[index_nk] == 0]
                lasso_pi <- lasso_var_selection(X[index_nk, ], A[index_nk], task = 1)
                selected_vars_pi <- lasso_pi$selected_vars
                pi_rf = caret::train(x = X[index_nk, selected_vars_pi, drop = FALSE], y = A_fac[index_nk], 
                  method = "ranger", trControl = ctrl_pi, tuneGrid = generate_tune_grid(length(index_nk), 
                    length(selected_vars_pi), 1), metric = "ROC", importance = "none", num.trees = generate_num_tree(length(index_nk)))
                prediction_pi = predict(pi_rf$finalModel, data = X[index, selected_vars_pi, drop = FALSE])
                pred_pi_a[4, index] = pmax(pmin(prediction_pi$predictions[, 2], 0.99), 0.01)
                pred_pi_b[4, index] = pred_pi_a[4, index]
                mtry_pi = c(mtry_pi, pi_rf$bestTune$mtry)
                minsize_pi = c(minsize_pi, pi_rf$bestTune$min.node.size)
                selected_M <- select_different_M_components_ks(X[index_nk, ], M[index_nk, ], A[index_nk], 
                  max_selected = sdr_max_M, min_selected = 1)
                x_vars <- multi_y_var_selection(X[index_nk, ], M[index_nk, selected_M, drop = FALSE])$union_vars
                basis_sel <- select_el_drm_basis(X = X[index_nk, x_vars, drop = FALSE], M = M[index_nk, 
                  selected_M, drop = FALSE], A = A[index_nk], stability_reps = 10, subsample_frac = 0.8, 
                  cv_folds = 5, alpha = 0.5, stability_cutoff = 0.8, seed = seed + 20000 * t + 200 * 
                    k)
                selected_basis <- basis_sel$selected_basis
                log_info(paste("QR2 basis fold", k, "; selected_M =", paste(selected_M, collapse = ","), 
                  "; selected_basis =", paste(selected_basis, collapse = ","), "; successful_reps =", 
                  basis_sel$successful_reps, "; quadratic_freq =", paste(names(basis_sel$quadratic_frequency), 
                    round(basis_sel$quadratic_frequency, 2), sep = ":", collapse = ",")))
                f_el_drm <- el_drm_fit_density_ratio(X[index_nk, x_vars, drop = FALSE], M[index_nk, selected_M, 
                  drop = FALSE], A[index_nk], basis_names = selected_basis, ridge_grid = sdr_ridge_grid, 
                  ridge_cv_folds = 5, ridge_cv_seed = seed + 10000 * t + 100 * k + 2)
                log_ratio_index <- el_drm_predict_log_ratio(f_el_drm, X[index, x_vars, drop = FALSE], 
                  M[index, selected_M, drop = FALSE])
                ref_log_q_sdr2_full <- f_el_drm$full_train_log_ratio[A[index_nk] == 1]
                pred_q[4, index] <- normalized_truncate_q(log_q = log_ratio_index, ref_log_q = ref_log_q_sdr2_full)
                log_info(paste("QR2 fold", k, "ridge_lambda =", f_el_drm$ridge_lambda, "x_dim =", f_el_drm$x_dim, 
                  "q_dim =", f_el_drm$q_dim, "theta_dim =", length(f_el_drm$theta)))
                lasso_mu <- lasso_var_selection(cbind(X[index_nk_A, ], M[index_nk_A, ]), Y[index_nk_A], 
                  task = 0)
                selected_vars_mu <- lasso_mu$selected_vars
                mu_rf = caret::train(x = cbind(X[index_nk_A, ], M[index_nk_A, ])[, selected_vars_mu, 
                  drop = FALSE], y = Y[index_nk_A], method = "ranger", tuneGrid = generate_tune_grid(length(index_nk_A), 
                  length(selected_vars_mu), 0), trControl = ctrl_mu, importance = "none", num.trees = generate_num_tree(length(index_nk_A)))
                pred_mu[4, index] = predict(mu_rf$finalModel, data = cbind(X[index, ], M[index, ])[, 
                  selected_vars_mu, drop = FALSE])$predictions
                mtry_mu = c(mtry_mu, mu_rf$bestTune$mtry)
                minsize_mu = c(minsize_mu, mu_rf$bestTune$min.node.size)
                prediction_mu_1 = predict(mu_rf$finalModel, data = cbind(X[index_nk_A0, ], M[index_nk_A0, 
                  ])[, selected_vars_mu, drop = FALSE])
                lasso_tau1 <- lasso_var_selection(X[index_nk_A0, ], prediction_mu_1$predictions, task = 0)
                selected_vars_tau1 <- lasso_tau1$selected_vars
                tau1_rf = caret::train(x = X[index_nk_A0, selected_vars_tau1, drop = FALSE], y = prediction_mu_1$predictions, 
                  method = "ranger", tuneGrid = generate_tune_grid(length(index_nk_A0), length(selected_vars_tau1), 
                    0), trControl = ctrl_mu, importance = "none", num.trees = generate_num_tree(length(index_nk_A0)))
                pred_tau1[4, index] = predict(tau1_rf$finalModel, data = X[index, selected_vars_tau1, 
                  drop = FALSE])$predictions
                mtry_tau1 = c(mtry_tau1, tau1_rf$bestTune$mtry)
                minsize_tau1 = c(minsize_tau1, tau1_rf$bestTune$min.node.size)
                prediction_q_1 <- normalized_truncate_q(log_q = ref_log_q_sdr2_full, ref_log_q = ref_log_q_sdr2_full)
                prediction_mu_2_preds <- predict_caret_regression_final(mu_rf, newdata = cbind(X[index_nk_A, 
                  ], M[index_nk_A, ])[, selected_vars_mu, drop = FALSE])
                prediction_tau1_1 = predict(tau1_rf$finalModel, data = X[index_nk_A, selected_vars_tau1, 
                  drop = FALSE])
                pred_pseudoY = prediction_q_1 * (Y[index_nk_A] - prediction_mu_2_preds) + prediction_tau1_1$predictions
                lasso_tau <- lasso_var_selection(X[index_nk_A, ], pred_pseudoY, task = 0)
                selected_vars_tau <- lasso_tau$selected_vars
                tau_rf = caret::train(x = X[index_nk_A, selected_vars_tau, drop = FALSE], y = pred_pseudoY, 
                  method = "ranger", tuneGrid = generate_tune_grid(length(index_nk_A), length(selected_vars_tau), 
                    0), trControl = ctrl_mu, importance = "none", num.trees = generate_num_tree(length(index_nk_A)))
                pred_tau[4, index] = predict(tau_rf$finalModel, data = X[index, selected_vars_tau, drop = FALSE])$predictions
                mtry_tau = c(mtry_tau, tau_rf$bestTune$mtry)
                minsize_tau = c(minsize_tau, tau_rf$bestTune$min.node.size)
                mu_rf = butcher(mu_rf)
                tau1_rf = butcher(tau1_rf)
                tau_rf = butcher(tau_rf)
                rm(mu_rf, tau1_rf, tau_rf, prediction_mu_1, prediction_mu_2_preds, prediction_tau1_1)
                gc()
            }
            log_info(paste("rep", t, "QR2 complete"))
            cat(format(Sys.time()), sprintf(" | rep %03d | QR2 complete\n", t))
        }
        if (method == "O-TR") {
            set.seed(seed)
            for (k in 1:Kmb) {
                index = (1:n) + (k - 1) * n
                index_nk = (1:N)[-index]
                index_nk1 = index_nk[1:floor(length(index_nk)/2)]
                index_nk2 = index_nk[-(1:floor(length(index_nk)/2))]
                index_nk1_A = index_nk1[A[index_nk1] == 1]
                index_nk2_A0 = index_nk2[A[index_nk2] == 0]
                index_nk_A = index_nk[A[index_nk] == 1]
                index_nk_A0 = index_nk[A[index_nk] == 0]
                n_nk_A = length(index_nk_A)
                n_nk_A0 = length(index_nk_A0)
                lasso_pi <- lasso_var_selection(X[index_nk, ], A[index_nk], task = 1)
                selected_vars_pi <- lasso_pi$selected_vars
                pi_rf_TR = caret::train(x = X[index_nk, selected_vars_pi, drop = FALSE], y = A_fac[index_nk], 
                  method = "ranger", trControl = ctrl_pi, tuneGrid = generate_tune_grid(length(index_nk), 
                    length(selected_vars_pi), 1), metric = "ROC", importance = "none", num.trees = generate_num_tree(length(index_nk)))
                prediction_pi_TR = predict(pi_rf_TR$finalModel, data = X[index, selected_vars_pi, drop = FALSE])
                pred_pi_a[6, index] = pmax(pmin(prediction_pi_TR$predictions[, 2], 0.99), 0.01)
                pred_pi_b[6, index] = pred_pi_a[6, index]
                pred_pi_a[7, index] = pred_pi_a[6, index]
                pred_pi_b[7, index] = pred_pi_a[6, index]
                pred_pi_a[9, index] = pred_pi_a[6, index]
                pred_pi_b[9, index] = pred_pi_a[6, index]
                mtry_pi = c(mtry_pi, pi_rf_TR$bestTune$mtry)
                minsize_pi = c(minsize_pi, pi_rf_TR$bestTune$min.node.size)
                pi_rf_TR = butcher(pi_rf_TR)
                rm(pi_rf_TR, prediction_pi_TR)
                gc()
                lasso_p <- lasso_var_selection(cbind(X[index_nk, ], M[index_nk, ]), A[index_nk], task = 1)
                selected_vars_p <- lasso_p$selected_vars
                p_rf_TR = caret::train(x = cbind(X[index_nk, ], M[index_nk, ])[, selected_vars_p, drop = FALSE], 
                  y = A_fac[index_nk], method = "ranger", trControl = ctrl_pi, tuneGrid = generate_tune_grid(length(index_nk), 
                    length(selected_vars_p), 1), metric = "ROC", importance = "none", num.trees = generate_num_tree(length(index_nk)))
                pred_p[4, index] = pmax(pmin(predict(p_rf_TR$finalModel, data = cbind(X[index, ], M[index, 
                  ])[, selected_vars_p, drop = FALSE])$predictions[, 2], 0.99), 0.01)
                pred_q[7, index] = (1 - pred_p[4, index])/pred_p[4, index] * pred_pi_a[7, index]/(1 - 
                  pred_pi_a[7, index])
                pred_p[6, index] = pred_p[4, index]
                pred_q[9, index] = pred_q[7, index]
                mtry_p = c(mtry_p, p_rf_TR$bestTune$mtry)
                minsize_p = c(minsize_p, p_rf_TR$bestTune$min.node.size)
                p_rf_TR = butcher(p_rf_TR)
                rm(p_rf_TR)
                lasso_mu <- lasso_var_selection(cbind(X[index_nk1_A, ], M[index_nk1_A, ]), Y[index_nk1_A], 
                  task = 0)
                selected_vars_mu <- lasso_mu$selected_vars
                mu_rf_TR = caret::train(x = cbind(X[index_nk1_A, ], M[index_nk1_A, ])[, selected_vars_mu, 
                  drop = FALSE], y = Y[index_nk1_A], method = "ranger", trControl = ctrl_mu, importance = "none", 
                  num.trees = generate_num_tree(length(index_nk1_A)), tuneGrid = generate_tune_grid(length(index_nk1_A), 
                    length(selected_vars_mu), 0))
                mtry_mu = c(mtry_mu, mu_rf_TR$bestTune$mtry)
                minsize_mu = c(minsize_mu, mu_rf_TR$bestTune$min.node.size)
                pred_mu[7, index] = predict(mu_rf_TR$finalModel, data = cbind(X[index, ], M[index, ])[, 
                  selected_vars_mu, drop = FALSE])$predictions
                prediction_mu_new1_TR = predict(mu_rf_TR$finalModel, data = cbind(X[index_nk2_A0, ], 
                  M[index_nk2_A0, ])[, selected_vars_mu, drop = FALSE])
                lasso_tau1 <- lasso_var_selection(X[index_nk2_A0, ], prediction_mu_new1_TR$predictions, 
                  task = 0)
                selected_vars_tau1 <- lasso_tau1$selected_vars
                tau1_rf_TR = caret::train(x = X[index_nk2_A0, selected_vars_tau1, drop = FALSE], y = prediction_mu_new1_TR$predictions, 
                  method = "ranger", trControl = ctrl_mu, tuneGrid = generate_tune_grid(length(index_nk2_A0), 
                    length(selected_vars_tau1), 0), importance = "none", num.trees = generate_num_tree(length(index_nk2_A0)))
                pred_tau1[7, index] = predict(tau1_rf_TR$finalModel, data = X[index, selected_vars_tau1, 
                  drop = FALSE])$predictions
                mtry_tau1 = c(mtry_tau1, tau1_rf_TR$bestTune$mtry)
                minsize_tau1 = c(minsize_tau1, tau1_rf_TR$bestTune$min.node.size)
            }
        }
        pred_tau[7, ] <- pred_tau1[7, ]
        ids <- c(3L, 4L, 7L)
        method_names <- c("QR1", "QR2", "O-TR", "Oracle")
        scores <- rbind(QR1 = A * pred_q[3, ]/pred_pi_a[3, ] * (Y - pred_mu[3, ]) + (1 - A)/(1 - pred_pi_b[3, 
            ]) * (pred_mu[3, ] - pred_tau1[3, ]) + A/pred_pi_a[3, ] * (pred_tau1[3, ] - pred_tau[3, ]) + 
            pred_tau[3, ], QR2 = A * pred_q[4, ]/pred_pi_a[4, ] * (Y - pred_mu[4, ]) + (1 - A)/(1 - pred_pi_b[4, 
            ]) * (pred_mu[4, ] - pred_tau1[4, ]) + A/pred_pi_a[4, ] * (pred_tau1[4, ] - pred_tau[4, ]) + 
            pred_tau[4, ], `O-TR` = A/(1 - pred_pi_a[7, ]) * (1 - pred_p[4, ])/pred_p[4, ] * (Y - pred_mu[7, 
            ]) + (1 - A)/(1 - pred_pi_a[7, ]) * (pred_mu[7, ] - pred_tau1[7, ]) + pred_tau1[7, ])
    }
    if (group == "qr") {
        keep <- match(method, method_names)
        ids <- ids[keep]
        method_names <- c(method, "Oracle")
        scores <- scores[method, , drop = FALSE]
    }
    oracle_score <- A * q/pi * (Y - mu) + (1 - A)/(1 - pi) * (mu - tau) + tau
    scores <- rbind(scores, Oracle = oracle_score)
    pred_theta <- rowMeans(scores)
    pred_sig2 <- rowMeans((scores - pred_theta)^2)
    stopifnot(all(is.finite(pred_theta)), all(is.finite(pred_sig2)), all(pred_sig2 >= 0))
    if (abs(tail(pred_theta, 1) - meta$oracle_theta[t]) > 1e-10 || abs(tail(pred_sig2, 1) - meta$oracle_variance[t]) > 
        1e-08) 
        stop("Oracle mismatch in rep ", t)
    take <- function(x) {
        z <- x[ids, , drop = FALSE]
        rownames(z) <- method_names[-length(method_names)]
        z
    }
    nuisance_mse <- t(vapply(ids, function(j) c(pi_a = mean((pred_pi_a[j, ] - pi)^2), pi_b = mean((pred_pi_b[j, 
        ] - pi)^2), q = mean((pred_q[j, A == 1] - q[A == 1])^2), mu = mean((pred_mu[j, ] - mu)^2), tau_n = mean((pred_tau1[j, 
        ] - tau)^2), tau = mean((pred_tau[j, ] - tau)^2)), numeric(6)))
    rownames(nuisance_mse) <- method_names[-length(method_names)]
    nuisance_mse <- rbind(nuisance_mse, Oracle = rep(0, 6))
    record <- list(replication_id = as.integer(t), case_id = meta$config$case_id, N = N, seed = seed, 
        theta = theta, method_names = method_names, pred_theta = pred_theta, pred_sig2 = pred_sig2, nuisance_mse = nuisance_mse, 
        data_md5 = meta$cache_manifest$md5[t])
    validate_record(record, meta, method, t)
    atomic_rds(record, file.path(run_dir, "replications", sprintf("rep_%03d.rds", t)))
    cat(format(Sys.time()), sprintf(" | rep %03d | %s complete\n", t, group))
    list(pred_theta = pred_theta, pred_sig2 = pred_sig2, nuisance_mse = nuisance_mse)
}

base_parameters <- function (config) 
{
    d <- config$d2
    c_M <- c(0.5, 0.5, rep(sqrt(0.5/(d - 2L)), d - 2L))
    gamma <- c(sqrt(0.87), -sqrt(0.87), rep(0, d - 2L))
    if (config$gamma_offset != 0) {
        gamma[1:2] <- gamma[1:2] + config$gamma_offset
        gamma <- gamma * sqrt(1.74/sum(gamma^2))
    }
    contrast <- c(1, -1, rep(0, d - 2L))/sqrt(2)
    rotation_vector <- numeric(d)
    if (config$axis_aligned_mu) {
        axis <- c(1, rep(0, d - 1L))
        rotation_vector <- gamma/sqrt(1.74) - axis
        rotation_vector <- rotation_vector/sqrt(sum(rotation_vector^2))
        c_M <- c_M - 2 * rotation_vector * sum(rotation_vector * c_M)
        contrast <- contrast - 2 * rotation_vector * sum(rotation_vector * contrast)
        gamma <- sqrt(1.74) * axis
    }
    beta_tau_X <- numeric(config$d1)
    beta_tau_X[100L] <- 0.4 * sqrt(1.74)
    vX <- 1 - 4 * dnorm(2)/(pnorm(2) - pnorm(-2))
    beta_pi_X <- numeric(config$d1)
    beta_pi_X[1] <- 0.8
    beta_pi_X[4:53] <- 0.7/sqrt(50 * vX)
    shared <- 1 - config$shifted_probability
    qcap <- if (shared > 0) 
        1/shared
    else Inf
    stopifnot(abs(sum(c_M^2) - 1) < 1e-12, abs(sum(gamma^2) - 1.74) < 1e-12, config$d1 >= 100L)
    list(configuration = config, c_M = c_M, delta = c_M, beta_mu_M = gamma, beta_mu_X = beta_tau_X, beta_tau_X = beta_tau_X, 
        B_scale = 0, beta_pi_X = beta_pi_X, beta_pi_i = 0.25, rho = 0, sigma = config$sigma_Y, trunc_variance = vX, 
        treated_scale = 1, rare_probability = 0.02, pi_bounds = c(0.1, 0.9), q_true_upper_bound = qcap, 
        q_second_moment_A1_upper_bound = if (is.finite(qcap)) qcap else exp((config$shift_mid + config$shift_span)^2), 
        mu_mediator_variance = 1.74, mu_mediator_support = as.integer(sum(gamma != 0)), tau_slope_support = 1L, 
        q_projection_support = d, pi_linear_support = 51L, shifted_probability = config$shifted_probability, 
        shared_probability = shared, outcome_noise = config$outcome_noise, outcome_noise_variance = config$sigma_Y^2, 
        mediator_contrast = contrast, mediator_rotation_vector = rotation_vector, axis_aligned_mu = config$axis_aligned_mu, 
        mediator_region_probability = (pnorm(2) - pnorm(1))/(pnorm(2) - pnorm(-2)), mediator_variance_floor = 0.2, 
        mediator_variance_mode = config$mediator_variance_mode, outcome_noise_fourth_moment = (if (config$outcome_noise == 
            "student_t5") 9 else 3) * config$sigma_Y^4, extra_X_start = 101L, true_mean_transport_overlap = sum(gamma * 
            c_M), true_pi_tau_shared_coordinate = FALSE)
}

weight_features <- function (X, config, parameters) 
{
    L <- rowSums(X[, 4:53, drop = FALSE])/sqrt(50 * parameters$trunc_variance)
    H <- 0.75 * X[, 2] * X[, 3] + 0.5 * (X[, 54]^2 - parameters$trunc_variance)
    list(L = L, H = H)
}

propensity_truth <- function (X, config, parameters) 
{
    f <- weight_features(X, config, parameters)
    0.1 + 0.8 * plogis(0.25 + 0.8 * X[, 1] + 0.7 * f$L + 0.75 * f$H + 0.3 * sin(f$L))
}

mediator_geometry <- function (X, config, parameters) 
{
    region <- switch(config$mediator_variance_mode, constant = rep(FALSE, nrow(X)), high_pi = X[, 1] > 
        1, low_pi = X[, 1] < -1, stop("Unknown mediator variance mode."))
    variance <- if (config$mediator_variance_mode == "constant") 
        rep(1, nrow(X))
    else parameters$mediator_variance_floor + (1 - parameters$mediator_variance_floor) * region/parameters$mediator_region_probability
    list(region = region, variance_contrast = variance)
}

true_log_ratio <- function (shift, projection, config) 
{
    lr <- shift * projection - 0.5 * shift^2
    if (config$shifted_probability == 1) 
        return(-lr)
    x <- log1p(-config$shifted_probability)
    y <- log(config$shifted_probability) + lr
    -(pmax(x, y) + log1p(exp(-abs(x - y))))
}

auxiliary_variance <- function (X, parameters) 
{
    parameters$auxiliary_variance_normalizer/(1 + parameters$auxiliary_kappa * X[, 1L])
}

required_packages <- function (methods) 
{
    unique(c("truncnorm", if (any(methods != "Oracle")) c("glmnet", "RCAL", "matrixStats", "logger"), 
        if (any(methods %in% c("QR1", "QR2", "TR", "O-TR"))) c("caret", "ranger", "randomForestSRC", 
            "butcher", "corpcor"), if (any(methods %in% c("QR1", "QR2"))) "MASS", if ("TR" %in% methods) "mvtnorm"))
}

load_fit_packages <- function (methods) 
{
    one_thread()
    pk <- required_packages(methods)
    missing <- pk[!vapply(pk, requireNamespace, logical(1), quietly = TRUE)]
    if (length(missing)) 
        stop("Missing packages: ", paste(missing, collapse = ", "))
    suppressPackageStartupMessages(invisible(lapply(pk, library, character.only = TRUE)))
    invisible(NULL)
}


# Data-generating process: Settings and coefficient vectors.
.simu7_config <- function(N) {
  stopifnot(length(N)==1L,is.finite(N),N>0,N==as.integer(N),N%%25L==0L)
  list(N=as.integer(N),seed=2026L,d1=100L,d2=50L,tmax=100L,Kmb=5L,case_id=49L,
    pi_correct=FALSE,q_correct=FALSE,simple_mean=TRUE,rare_X=FALSE,
    shift_mid=1.2,shift_span=.3,dense_q=FALSE,control_mean_X1=0,
    noise_mode='gaussian_measurement',mean_sign=-1,tau_coordinate=100L,
    shifted_probability=1,gamma_offset=.1,outcome_noise='normal',sigma_Y=.2,
    mediator_variance_mode='high_pi',axis_aligned_mu=TRUE,
    auxiliary_dimension=20L,auxiliary_kappa=.45)
}
.simu7_parameters <- function(config) {
  p<-base_parameters(config)
  c1<-.35
  p$c_M<-c(c1,rep(sqrt((1-c1^2)/2),2L),rep(0,config$d2-3L))
  p$delta<-p$c_M
  p$mediator_contrast<-(c(1,rep(0,config$d2-1L))-c1*p$c_M)/sqrt(1-c1^2)
  p$q_projection_support<-3L
  p$true_mean_transport_overlap<-sum(p$beta_mu_M*p$c_M)
  p$auxiliary_indices<-4:23
  p$auxiliary_kappa<-config$auxiliary_kappa
  k<-p$auxiliary_kappa
  p$auxiliary_variance_normalizer<-(pnorm(2)-pnorm(-2))/integrate(
    function(z)dnorm(z)/(1+k*z),-2,2,rel.tol=1e-12)$value
  p$auxiliary_variance_bounds<-p$auxiliary_variance_normalizer/c(1+2*k,1-2*k)
  stopifnot(all(p$c_M[p$auxiliary_indices]==0),
    all(p$mediator_contrast[p$auxiliary_indices]==0),all(p$beta_mu_M[p$auxiliary_indices]==0))
  p
}
shift_truth <- function(X,config,parameters) {
  config$shift_mid+config$shift_span*tanh(.75*X[,2]+.5*X[,3])
}
.simu7_generate <- function(config,parameters,t) {
  data_seed<-as.integer(config$seed+100003L*t)
  auxiliary_seed<-as.integer(data_seed+900001L)
  RNGkind('Mersenne-Twister','Inversion','Rejection');set.seed(data_seed)
  N<-config$N;d2<-config$d2
  X<-matrix(truncnorm::rtruncnorm(N*100L,a=-2,b=2),N,100L)
  pi<-propensity_truth(X,config,parameters)
  A<-rbinom(N,1,pi);epsilon<-rnorm(N,sd=parameters$sigma)
  M<-matrix(rnorm(N*d2),N,d2)
  M<-M-2*tcrossprod(as.numeric(M%*%parameters$mediator_rotation_vector),parameters$mediator_rotation_vector)
  geo<-mediator_geometry(X,config,parameters)
  contrast_noise<-as.numeric(M%*%parameters$mediator_contrast)
  M<-M+tcrossprod((sqrt(geo$variance_contrast)-1)*contrast_noise,parameters$mediator_contrast)
  shifted<-rep(1L,N);a<-shift_truth(X,config,parameters)
  M<-M+tcrossprod(A*shifted*a,parameters$c_M)
  w<-as.numeric(M%*%parameters$c_M)
  log_q<-true_log_ratio(a,w,config);q<-exp(log_q)
  p<-plogis(qlogis(pi)-log_q)
  tau<-1+as.numeric(X%*%parameters$beta_tau_X)
  mu<-tau+as.numeric(M%*%parameters$beta_mu_M)
  d<-list(replication_id=as.integer(t),N=N,seed=config$seed,data_seed=data_seed,case_id=config$case_id,
    X=X,M=M,A=A,Y=mu-(1-A)+epsilon,pi=pi,p=p,q=q,mu=mu,tau=tau,tau1=tau,
    mediator_variance=geo$variance_contrast,mediator_region=geo$region,
    outcome_variance=rep(parameters$sigma^2,N),noise_region=rep(FALSE,N),shift=a,
    treated_scale=1,auxiliary_seed=auxiliary_seed,shifted_component=shifted)
  jj<-parameters$auxiliary_indices
  d$M[,jj]<-d$M[,jj,drop=FALSE]*sqrt(auxiliary_variance(d$X,parameters))
  stopifnot(all(is.finite(d$q)),all(d$q>0),all(d$pi>=.1&d$pi<=.9))
  d
}
# Construct the same complete finite pool for a full run or any requested subset.
# Only requested replications are stored on disk, never exported to every worker.
.simu7_prepare_dgp <- function(N, replication_ids=integer(), cache_dir=NULL) {
  config<-.simu7_config(N);parameters<-.simu7_parameters(config)
  tau_all<-numeric(N*100L);oracle_theta<-oracle_variance<-numeric(100L)
  manifest<-data.frame(replication_id=1:100,md5=rep('',100L))
  if(length(replication_ids)&&is.null(cache_dir))stop('A cache directory is required for requested data.')
  if(!is.null(cache_dir))dir.create(cache_dir,recursive=TRUE,showWarnings=FALSE)
  for(id in 1:100){
    d<-.simu7_generate(config,parameters,id)
    tau_all[(id-1L)*N+seq_len(N)]<-d$tau
    score<-with(d,A*q/pi*(Y-mu)+(1-A)/(1-pi)*(mu-tau)+tau)
    oracle_theta[id]<-mean(score);oracle_variance[id]<-mean((score-mean(score))^2)
    if(id%in%replication_ids){
      path<-file.path(cache_dir,sprintf('rep_%03d.rds',id));atomic_rds(d,path,compress=FALSE)
      manifest$md5[id]<-unname(tools::md5sum(path))
    }
  }
  list(config=config,parameters=parameters,theta=mean(tau_all),oracle_theta=oracle_theta,
    oracle_variance=oracle_variance,cache_manifest=manifest,grids=tuning_grids())
}




# Execute one replication. Every method retains its original independent seed reset.
.simu7_run_replication <- function(id, slot, meta, cache_dir, run_dir) {
  load_fit_packages(.simu7_method_keys)
  input<-file.path(cache_dir,sprintf('rep_%03d.rds',id))
  if(!identical(unname(tools::md5sum(input)),meta$cache_manifest$md5[id]))stop('Data checksum mismatch.')
  d<-readRDS(input)
  estimates<-variances<-setNames(numeric(8),.simu7_method_names)
  estimates['Oracle']<-meta$oracle_theta[id];variances['Oracle']<-meta$oracle_variance[id]
  temporary<-file.path(run_dir,'.working',sprintf('rep_%03d',id))
  dir.create(temporary,recursive=TRUE,showWarnings=FALSE)
  for(key in setdiff(.simu7_method_keys,'Oracle')){
    method<-method_label(key);md<-file.path(temporary,key);rd<-file.path(md,'replications')
    dir.create(rd,recursive=TRUE,showWarnings=FALSE)
    logger::log_appender(logger::appender_file(meta$log_file))
    logger::log_info(sprintf('Replication %03d/100 | %s | started',id,method))
    if(key%in%c('MQR1','MQR2'))fit_parametric_method(d,meta,key,meta$log_file,rd)
    else fit_comparator_method(id,meta,cache_dir,md,key)
    record<-readRDS(file.path(rd,sprintf('rep_%03d.rds',id)))
    validate_record(record,meta,key,id)
    estimates[method]<-record$pred_theta[method];variances[method]<-record$pred_sig2[method]
    logger::log_info(sprintf('Replication %03d/100 | %s | complete',id,method))
    rm(record);invisible(gc(FALSE))
  }
  if(any(!is.finite(estimates))||any(!is.finite(variances))||any(variances<0))stop('Invalid method result.')
  record<-list(source_id=paste0('simu_7_N',meta$config$N),replication_id=as.integer(id),
    N=meta$config$N,seed=meta$config$seed,theta=meta$theta,method_names=.simu7_method_names,
    pred_theta=estimates,pred_sig2=variances,data_md5=meta$cache_manifest$md5[id])
  atomic_rds(record,file.path(run_dir,'replications',sprintf('rep_%03d.rds',id)))
  unlink(temporary,recursive=TRUE)
  unlink(input)
  logger::log_info(sprintf('Replication %03d/100 | results | saved',id))
  invisible(id)
}

# Public interface shared with simu_1 through simu_6.
run_experiment <- function(rerun_config,rerun_replication_ids,rerun_workers,rerun_output_dir) {
  if(length(rerun_config)!=1L||!rerun_config%in%c(1000,1800,2500,5000))stop('Unsupported N.')
  if(!length(rerun_replication_ids)||anyNA(rerun_replication_ids)||
    any(rerun_replication_ids!=as.integer(rerun_replication_ids))||
    any(!rerun_replication_ids%in%1:100)||anyDuplicated(rerun_replication_ids))stop('Use distinct replications from1:100.')
  if(length(rerun_workers)!=1L||!is.finite(rerun_workers)||rerun_workers<1||
    rerun_workers!=as.integer(rerun_workers))stop('workers must be a positive integer.')
  ids<-as.integer(rerun_replication_ids)
  layout<-.prepare_run_output('simu_7',paste0('n_',rerun_config),ids,rerun_workers,rerun_output_dir)
  run_dir<-layout$run_directory;log_file<-layout$log_file
  old_rng_kind<-RNGkind()
  old_rng<-if(exists('.Random.seed',.GlobalEnv,inherits=FALSE))get('.Random.seed',.GlobalEnv)else NULL
  on.exit({do.call(RNGkind,as.list(old_rng_kind));if(is.null(old_rng)){
    if(exists('.Random.seed',.GlobalEnv,inherits=FALSE))rm('.Random.seed',envir=.GlobalEnv)
  }else assign('.Random.seed',old_rng,.GlobalEnv)},add=TRUE)
  load_fit_packages(.simu7_method_keys)
  cache_dir<-file.path(run_dir,'.working','data_cache')
  cat('Preparing the finite-pool target and disk-backed replication data.\n');flush.console()
  meta<-.simu7_prepare_dgp(as.integer(rerun_config),ids,cache_dir)
  meta$log_file<-log_file
  cat('Finite-pool theta: ',format(meta$theta,digits=17),'\nProgress log: ',log_file,'\n',sep='');flush.console()
  cl<-parallel::makePSOCKcluster(min(as.integer(rerun_workers),length(ids)),outfile=log_file)
  on.exit(if(!is.null(cl))try(parallel::stopCluster(cl),silent=TRUE),add=TRUE)
  env<-environment(run_experiment);fn<-ls(env,all.names=TRUE)
  fn<-fn[vapply(fn,function(n)is.function(get(n,env)),FALSE)]
  parallel::clusterExport(cl,c(fn,'.simu7_method_names','.simu7_method_keys'),envir=env)
  plan<-worker_plan(ids,length(cl))
  schedule_fixed_reps(cl,plan,.simu7_run_replication,args=list(meta=meta,cache_dir=cache_dir,run_dir=run_dir))
  parallel::stopCluster(cl);cl<-NULL
  rows<-lapply(ids,function(id)readRDS(file.path(layout$replication_directory,sprintf('rep_%03d.rds',id))))
  pred_theta<-do.call(rbind,lapply(rows,`[[`,'pred_theta'))
  pred_sig2<-do.call(rbind,lapply(rows,`[[`,'pred_sig2'))
  dimnames(pred_theta)<-dimnames(pred_sig2)<-list(as.character(ids),.simu7_method_names)
  theta<-meta$theta;N<-meta$config$N
  error<-pred_theta-theta;se<-sqrt(pred_sig2/N)
  bias<-apply(error,2,median);RMSE<-sqrt(colMeans(error^2))
  AL<-apply(2*qnorm(.975)*se,2,median);AC<-colMeans(abs(error)<=qnorm(.975)*se)
  ESD<-apply(pred_theta,2,function(x)1.4826*median(abs(x-median(x))));ASD<-apply(se,2,median)
  sig2<-apply(pred_sig2,2,median)
  paper_method_names<-c('Oracle','TR','O-TR','QR1','QR2',"O-TR'",'MQR1','MQR2')
  table_summary<-data.frame(method=.simu7_method_names,bias=unname(bias),RMSE=unname(RMSE),
    AL=unname(AL),AC=unname(AC),ESD=unname(ESD),ASD=unname(ASD))
  table_summary<-table_summary[match(paper_method_names,table_summary$method),];rownames(table_summary)<-NULL
  manuscript_table<-table_summary[,c('method','bias','RMSE','AL','AC')]
  for(k in c('bias','RMSE','AL','AC'))manuscript_table[[k]]<-round(manuscript_table[[k]],if(k=='AC')2 else 3)
  result<-list(schema_version='1.1',source_id=paste0('simu_7_N',N),group=7L,config=as.numeric(N),
    source_path=paste0('configurations/simu_7_N',N,'.R'),replication_ids=ids,tmax=100L,N=N,
    theta=theta,theta_definition='mean(tau_all)',seed=meta$config$seed,method_names=.simu7_method_names,
    method_columns=seq_len(8L),paper_method_names=paper_method_names,num_totest=8L,
    pred_theta=pred_theta,pred_sig2=pred_sig2,results_par=cbind(pred_theta,pred_sig2),
    sig2=sig2,bias=bias,RMSE=RMSE,AL=AL,AC=AC,ESD=ESD,ASD=ASD,
    table_summary=table_summary,manuscript_table=manuscript_table,run_directory=run_dir,
    log_file=log_file,session_info=sessionInfo())
  colnames(result$results_par)<-c(paste0('pred_theta:',.simu7_method_names),paste0('pred_sig2:',.simu7_method_names))
  path<-file.path(layout$results_directory,'results.RData')
  save(list=names(result),envir=list2env(result,parent=emptyenv()),file=path)
  write.csv(table_summary,file.path(layout$results_directory,'summary.csv'),row.names=FALSE)
  unlink(file.path(run_dir,'.working'),recursive=TRUE)
  print(table_summary,row.names=FALSE);message('Final result: ',path)
  invisible(path)
}

rerun_result_files <- setNames(vector('list',length(N)),as.character(N))
if(!isTRUE(getOption('qr.archive.no_autorun',FALSE))){
  for(.rerun_configuration in N){
    rerun_result_files[[as.character(.rerun_configuration)]]<-run_experiment(
      .rerun_configuration,replication_ids,workers,output_dir)
    gc(verbose=FALSE)
  }
}


