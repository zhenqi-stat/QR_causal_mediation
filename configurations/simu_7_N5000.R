# Setting(vii): data-generating process only; run estimation via simu_7.R.
N <- 5000L
replication_ids <- 1:100

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

.simu7_config <- function (N) 
{
    stopifnot(length(N) == 1L, is.finite(N), N > 0, N == as.integer(N), N%%25L == 0L)
    list(N = as.integer(N), seed = 2026L, d1 = 100L, d2 = 50L, tmax = 100L, Kmb = 5L, case_id = 49L, 
        pi_correct = FALSE, q_correct = FALSE, simple_mean = TRUE, rare_X = FALSE, shift_mid = 1.2, shift_span = 0.3, 
        dense_q = FALSE, control_mean_X1 = 0, noise_mode = "gaussian_measurement", mean_sign = -1, tau_coordinate = 100L, 
        shifted_probability = 1, gamma_offset = 0.1, outcome_noise = "normal", sigma_Y = 0.2, mediator_variance_mode = "high_pi", 
        axis_aligned_mu = TRUE, auxiliary_dimension = 20L, auxiliary_kappa = 0.45)
}

.simu7_parameters <- function (config) 
{
    p <- base_parameters(config)
    c1 <- 0.35
    p$c_M <- c(c1, rep(sqrt((1 - c1^2)/2), 2L), rep(0, config$d2 - 3L))
    p$delta <- p$c_M
    p$mediator_contrast <- (c(1, rep(0, config$d2 - 1L)) - c1 * p$c_M)/sqrt(1 - c1^2)
    p$q_projection_support <- 3L
    p$true_mean_transport_overlap <- sum(p$beta_mu_M * p$c_M)
    p$auxiliary_indices <- 4:23
    p$auxiliary_kappa <- config$auxiliary_kappa
    k <- p$auxiliary_kappa
    p$auxiliary_variance_normalizer <- (pnorm(2) - pnorm(-2))/integrate(function(z) dnorm(z)/(1 + k * 
        z), -2, 2, rel.tol = 1e-12)$value
    p$auxiliary_variance_bounds <- p$auxiliary_variance_normalizer/c(1 + 2 * k, 1 - 2 * k)
    stopifnot(all(p$c_M[p$auxiliary_indices] == 0), all(p$mediator_contrast[p$auxiliary_indices] == 0), 
        all(p$beta_mu_M[p$auxiliary_indices] == 0))
    p
}

shift_truth <- function(X, config, parameters) {
    config$shift_mid + config$shift_span * tanh(0.75 * X[, 2] + 0.5 * X[, 3])
}

.simu7_generate <- function (config, parameters, t) 
{
    data_seed <- as.integer(config$seed + 100003L * t)
    auxiliary_seed <- as.integer(data_seed + 900001L)
    RNGkind("Mersenne-Twister", "Inversion", "Rejection")
    set.seed(data_seed)
    N <- config$N
    d2 <- config$d2
    X <- matrix(truncnorm::rtruncnorm(N * 100L, a = -2, b = 2), N, 100L)
    pi <- propensity_truth(X, config, parameters)
    A <- rbinom(N, 1, pi)
    epsilon <- rnorm(N, sd = parameters$sigma)
    M <- matrix(rnorm(N * d2), N, d2)
    M <- M - 2 * tcrossprod(as.numeric(M %*% parameters$mediator_rotation_vector), parameters$mediator_rotation_vector)
    geo <- mediator_geometry(X, config, parameters)
    contrast_noise <- as.numeric(M %*% parameters$mediator_contrast)
    M <- M + tcrossprod((sqrt(geo$variance_contrast) - 1) * contrast_noise, parameters$mediator_contrast)
    shifted <- rep(1L, N)
    a <- shift_truth(X, config, parameters)
    M <- M + tcrossprod(A * shifted * a, parameters$c_M)
    w <- as.numeric(M %*% parameters$c_M)
    log_q <- true_log_ratio(a, w, config)
    q <- exp(log_q)
    p <- plogis(qlogis(pi) - log_q)
    tau <- 1 + as.numeric(X %*% parameters$beta_tau_X)
    mu <- tau + as.numeric(M %*% parameters$beta_mu_M)
    d <- list(replication_id = as.integer(t), N = N, seed = config$seed, data_seed = data_seed, case_id = config$case_id, 
        X = X, M = M, A = A, Y = mu - (1 - A) + epsilon, pi = pi, p = p, q = q, mu = mu, tau = tau, tau1 = tau, 
        mediator_variance = geo$variance_contrast, mediator_region = geo$region, outcome_variance = rep(parameters$sigma^2, 
            N), noise_region = rep(FALSE, N), shift = a, treated_scale = 1, auxiliary_seed = auxiliary_seed, 
        shifted_component = shifted)
    jj <- parameters$auxiliary_indices
    d$M[, jj] <- d$M[, jj, drop = FALSE] * sqrt(auxiliary_variance(d$X, parameters))
    stopifnot(all(is.finite(d$q)), all(d$q > 0), all(d$pi >= 0.1 & d$pi <= 0.9))
    d
}

