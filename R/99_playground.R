StanModel <- R6::R6Class(
  classname = "StanModel",
  public = list(
    data = list(),
    params = list(),
    initialize = function(data, params) {
      self$data = data
      self$params = params
    },
    
    set_data = function() {
    },
    compile = function() {
      paste(
        private$compile_data(), 
        private$compile_params(),
        sep = "\n"
      )
    }
  ),
  active = list(
    stan_code = function() {
    }
  ),
  private = list(
    compile_data = function() {
      str <- vapply(self$data, write_data, character(1))
      paste0(
        "data {\n",
        paste0(paste0("  ", str, collapse = ";\n"), ";"),
        "\n}"
      )
    },
    compile_params = function() {
      str <- vapply(self$params, write_data, character(1))
      paste0(
        "parameters {\n",
        paste0(paste0("  ", str, collapse = ";\n"), ";"),
        "\n}"
      )
    },
    compile_transformed_parameters = function() {
    },
    
    compile_model = function() {
    }
  )
)


write_data <- function(d) {
  lower <- "lower" %in% names(d)
  upper <- "upper" %in% names(d)
  
  if (lower && upper) {
    bounds <- glue::glue("<lower={d$lower}, upper={d$upper}>")
  } else if (lower) {
    bounds <- glue::glue("<lower={d$lower}>")
  } else if (upper) {
    bounds <- glue::glue("<upper={d$upper}>")
  } else {
    bounds <- ""
  }
  
  if (d$type %in% c("int", "real")) {
    str <- glue::glue("{d$type}{bounds} {d$name}")
  } else if (d$type == "vector") {
    str <- glue::glue("vector[{d$rows}]{bounds} {d$name}")
  } else if (d$type == "matrix") {
    str <- glue::glue("matrix[{d$rows}, {d$cols}]{bounds} {d$name}")
  } else {
    stop("Unrecognized data type ", d$type)
  }
  str
}

data <- list(
  list(name = "n", type = "int"), 
  list(name = "p", type = "int"), 
  list(name = "g", type = "real"),
  list(name = "y", type = "vector", rows = "n"),
  list(name = "mu_b", type = "vector", rows = "p"),
  list(name = "mu_Sigma", type = "matrix", rows = "p", cols = "p"),
  list(name = "X", type = "matrix", "rows" = "n", "cols" = "p")
)

params <- list(
  list(name = "log_sigma_sq", type = "real"),
  list(name = "beta", type = "vector", rows = "p")
)

m <- StanModel$new(data, params)
cat(m$compile())

model_str <- "data {
  int n; 
  int p; 
  real g;
  vector[n] y; 
  matrix[n, p] X;
  matrix[p, p] Sigma; // (X'X)^(-1)
  vector[p] mu;
}

parameters {
  vector[p] beta;
  real log_sigma_sq;
}

transformed parameters {
  real<lower=0> sigma;
  sigma = sqrt(exp(log_sigma_sq));
}

model {
  log_sigma_sq ~ uniform(-5, 5);
  beta ~ multi_normal(mu, g * pow(sigma, 2) * Sigma);
  y ~ normal(to_vector(X * beta), sigma);
}
"

# n: number of rows
# p: numbver of variables
# g: g factor
# y: response variable
# mu_b: prior mean
# mu_Sigma: unadjusted prior covariance matrix
# X: design matrix