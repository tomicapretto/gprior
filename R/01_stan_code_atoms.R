library(glue)

# Add 2 spaces to the beginning of each line in a multi-line string
spacify = function(string) {
  string = unlist(strsplit(string, "\n", fixed=TRUE))
  paste0("  ", string, collapse = "\n")
}

BAD_NAMES = c(
  "for", "in", "while", "repeat", "until", "if", "then", "else", "true", 
  "false", "target", "int", "real", "vector", "simplex", "unit_vector", 
  "ordered", "positive_ordered", "row_vector", "matrix", "cholesky_factor_corr", 
  "cholesky_factor_cov", "corr_matrix", "cov_matrix", "functions", "model", 
  "data", "parameters", "quantities", "transformed", "generated",
  "var", "fvar", "STAN_MAJOR", "STAN_MINOR", "STAN_PATCH",
  "STAN_MATH_MAJOR", "STAN_MATH_MINOR", "STAN_MATH_PATCH",
  "alignas", "alignof", "and", "and_eq", "asm", "auto", "bitand", "bitor", 
  "bool", "break", "case", "catch", "char", "char16_t", "char32_t", "class", 
  "compl", "const", "constexpr", "const_cast", "continue", "decltype", 
  "default", "delete", "do", "double", "dynamic_cast", "else", "enum", 
  "explicit", "export", "extern", "false", "float", "for", "friend", "goto", 
  "if", "inline", "int", "long", "mutable", "namespace", "new", "noexcept",
  "not", "not_eq", "nullptr", "operator", "or", "or_eq", "private",
  "protected", "public", "register", "reinterpret_cast", "return",
  "short", "signed", "sizeof", "static", "static_assert", "static_cast",
  "struct", "switch", "template", "this", "thread_local", "throw", "true",
  "try", "typedef", "typeid", "typename", "union", "unsigned", "using",
  "virtual", "void", "volatile", "wchar_t", "while", "xor", "xor_eq"
)

Atom = R6::R6Class(
  classname = "Atom",
  public = list(
    name = NULL,
    lower = NULL,
    upper = NULL,
    initialize = function(name, lower = NULL, upper = NULL) {
      if (name %in% BAD_NAMES) stop(name, " is not allowed.")
      self$name = name
      self$lower = lower
      self$upper = upper
    },
    write_bounds = function() {
      if (!(is.null(self$lower) || is.null(self$upper))) {
        bounds = glue("<lower={self$lower}, upper={self$upper}>")
      } else if (!is.null(self$lower)) {
        bounds = glue("<lower={self$lower}>")
      } else if (!is.null(self$upper)) {
        bounds = glue("<upper={self$upper}>")
      } else {
        bounds = ""
      }
      bounds
    }
  )
)

Number = R6::R6Class(
  classname = "Number",
  inherit = Atom,
  public = list(
    type = NULL,
    initialize = function(name, type, ...) {
      super$initialize(name, ...)
      self$type = type
    },
    write = function() {
      bounds = self$write_bounds()
      glue("{self$type}{bounds} {self$name}")
    }
  )
)

Vector = R6::R6Class(
  classname = "Vector",
  inherit = Atom,
  public = list(
    rows = NULL,
    initialize = function(name, rows, ...) {
      super$initialize(name, ...)
      self$rows = rows
    },
    write = function() {
      bounds = self$write_bounds()
      glue("vector[{self$rows}]{bounds} {self$name}")
    }
  )
)

Matrix = R6::R6Class(
  classname = "Matrix",
  inherit = Atom,
  public = list(
    rows = NULL,
    cols = NULL,
    initialize = function(name, rows, cols, ...) {
      super$initialize(name, ...)
      self$rows = rows
      self$cols = cols
    },
    write = function() {
      bounds = self$write_bounds()
      glue("matrix[{self$rows}, {self$cols}]{bounds} {self$name}")
    }
  )
)

# Arrays are the only way to store sequences of integers, and some functions in 
# Stan, such as discrete distributions, require integer arguments.
Array = R6::R6Class(
  classname = "Array",
  inherit = Atom,
  public = list(
    type = NULL,
    shape = NULL,
    initialize = function(name, type, shape, ...) {
      super$initialize(name, ...)
      self$type = type
      self$shape = shape
    },
    write = function() {
      bounds = self$write_bounds()
      dims = paste0(self$shape, collapse = ", ")
      glue("array[{dims}] {self$type}{bounds} {self$name}")
    }
  )
)

SamplingStatement = R6::R6Class(
  classname = "SamplingStatement",
  public = list(
    name = NULL,
    expression = NULL,
    initialize = function(name, expression) {
      self$name = name
      self$expression = expression
    },
    write = function() {
      glue("{self$name} ~ {self$expression}")
    }
  )
)

Assignment = R6::R6Class(
  classname = "Assignment",
  public = list(
    name = NULL,
    expression = NULL,
    initialize = function(name, expression) {
      self$name = name
      self$expression = expression
    },
    write = function() {
      if (inherits(self$name, "Atom")) {
        lhs = self$name$write()
      } else {
        lhs = self$name
      }
      glue("{lhs} = {self$expression}")
    }
  )
)

ForLoop = R6::R6Class(
  classname = "ForLoop",
  public = list(
    var = NULL,
    seq = NULL,
    expressions = list(),
    initialize = function(var, seq, ...) {
      self$var = var
      self$seq = seq
      self$expressions = list(...)
    },
    write = function() {
      # Note this does not consider nested loops!
      expressions = lapply(
        self$expressions, 
        function(expr) glue("{expr$write()};")
      )
      expressions = paste0(expressions, collapse = "\n")
      glue(
        "for ({self$var} in {self$seq}) {{\n{spacify(expressions)}\n}}"
      )
    }
  )
)

Block = R6::R6Class(
  classname = "Block",
  public = list(
    name = NULL,
    expressions = list(),
    initialize = function(name, ...) {
      self$name = name
      self$expressions = list(...)
      
    }, 
    write = function() {
      expressions = lapply(
        self$expressions, function(expr) {
          if (inherits(expr, "ForLoop")) {
            glue("{expr$write()}")
          } else {
            glue("{expr$write()};")
          }
        })
      expressions = paste0(expressions, collapse = "\n")
      glue(
        "{self$name} {{\n{spacify(expressions)}\n}}"
      )
    }
  )
)

Data = R6::R6Class(
  classname = "Data",
  inherit = Block,
  public = list(
    initialize = function(...) {
      super$initialize("data", ...)
    }
  )
)

Parameters = R6::R6Class(
  classname = "Parameters",
  inherit = Block,
  public = list(
    initialize = function(...) {
      super$initialize("parameters", ...)
    }
  )
)

TransformedParameters = R6::R6Class(
  classname = "TransformedParameters",
  inherit = Block,
  public = list(
    initialize = function(...) {
      super$initialize("transformed parameters", ...)
    }
  )
)

Model = R6::R6Class(
  classname = "Model",
  inherit = Block,
  public = list(
    initialize = function(...) {
      super$initialize("model", ...)
    }
  )
)

GeneratedQuantities = R6::R6Class(
  classname = "GeneratedQuantities",
  inherit = Block,
  public = list(
    initialize = function(...) {
      super$initialize("generated quantities", ...)
    }
  )
) 

StanProgram = R6::R6Class(
  classname = "StanProgram",
  public = list(
    blocks = list(),
    initialize = function(...) {
      self$blocks = list(...)
    },
    write = function() {
      blocks = lapply(self$blocks, function(block) block$write())
      paste0(blocks, collapse = "\n")
    }
  )
)

# Examples ----------------------------------------------------------------
# number = Number$new("n", "int", 1, 2)
# number$write()
# 
# vector = Vector$new("master", 10, 0, 100)
# vector$write()
# 
# matrix = Matrix$new("X", 10, 3, 0, 100)
# matrix$write()
# 
# array = Array$new("arr", "int", c(20, 5), 0)
# array$write()
# 
# array = Array$new("arr", "real", c(20, 5), 0)
# array$write()
# 
# sampling_statement = SamplingStatement$new("y", "normal(mu, sigma)")
# sampling_statement$write()
# 
# assignment = Assignment$new("log_lik[i]", "normal_lpdf(y[i] | X[i] * beta, sigma)")
# assignment$write()
# 
# assignment = Assignment$new(
#   Number$new("sigma_p", "real"), 
#   "exponential_rng(1 / sd(y))"
# )
# assignment$write()
# 
# loop = ForLoop$new(
#   "i", "1:N",
#   Assignment$new("log_lik[i]", "normal_lpdf(y[i] | X[i] * beta, sigma)")
# )
# loop$write()
# 
# 
# program = StanProgram$new(
#   Data$new(
#     Number$new("n", "int"),
#     Number$new("p", "int"),
#     Number$new("g", "real"),
#     Vector$new("y", "n"),
#     Vector$new("mu_b", "p"),
#     Matrix$new("mu_Sigma", "p", "p"),
#     Matrix$new("X", "n", "p")
#   ),
#   Parameters$new(
#     Vector$new("beta", "p"),
#     Number$new("log_sigma_sq", "real")
#   ),
#   TransformedParameters$new(
#     Assignment$new(
#       Number$new("sigma", "real", lower=0),
#       "sqrt(exp(log_sigma_sq))"
#     )
#   ),
#   Model$new(
#     SamplingStatement$new(
#       "beta",
#       "multi_normal(mu_b, g * pow(sigma, 2) * mu_Sigma)"
#     ),
#     SamplingStatement$new(
#       "y",
#       "normal(to_vector(X * beta), sigma)"
#     )
#   ),
#   GeneratedQuantities$new(
#     Vector$new("log_lik", "n"),
#     ForLoop$new(
#       "i", "1:n",
#       Assignment$new("log_lik[i]", "normal_lpdf(y[i] | X[i] * beta, sigma)")
#     )
#   )
# )
# cat(program$write())