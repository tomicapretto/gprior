library(glue)

StanModelCode = R6::R6Class(
  classname = "StanModel",
  public = list(
    model = list(),
    initialize = function(model = list()) {
      self$model = model
    },
    add_line = function(block, line) {
      check_block_order = !block %in% names(self$model)
      name = line$name
      line$name = NULL
      self$model[[block]][[name]] = line
      if (check_block_order) {
        self$sort_blocks()
      }
    },
    add_lines = function(block, lines) {
      for (line in lines) {
        self$add_line(block, line)
      }
    },
    remove_line = function(block, name) {
      self$model[[block]][[name]] = NULL
    },
    remove_lines = function(block, name) {
      for (line in lines) {
        self$remove_line(block, line)
      }
    },
    remove_block = function(block) {
      self$model[[block]] = NULL
    },
    sort_blocks = function() {
      new_pos = match(private$BLOCK_ORDER, names(self$model))
      new_pos = new_pos[!is.na(new_pos)]
      self$model = self$model[new_pos]
    },
    sort_lines = function() {
      
      
    },
    make_stan_code = function() {
      code_list = mapply(
        function(block_name, block_components) {
          str = mapply(
            function(name, props) private$make_line(name, props),
            names(block_components),
            block_components,
            SIMPLIFY = FALSE
          )
          str = as.character(str)
          paste0(
            block_name, " {\n",
            paste0(add_space(split_lines(str)), collapse = "\n"),
            "\n}"
          )
        }, 
        names(self$model),
        self$model,
        SIMPLIFY = FALSE
      )
      paste(as.character(code_list), collapse = "\n")
    }
  ),
  private = list(
    
    BLOCK_ORDER = c(
      "data", "parameters", "transformed parameters", "model", 
      "generated quantities"
    ),
    
    make_line = function(name, props) {
      lower = "lower" %in% names(props)
      upper = "upper" %in% names(props)
      
      if (lower && upper) {
        bounds = glue::glue("<lower={props$lower}, upper={props$upper}>")
      } else if (lower) {
        bounds = glue::glue("<lower={props$lower}>")
      } else if (upper) {
        bounds = glue::glue("<upper={props$upper}>")
      } else {
        bounds = ""
      }
      
      if (props$type %in% c("int", "real")) {
        if ("dims" %in% names(props)) {
          str = glue::glue("{props$type}{bounds} {name}{props$dims}")
        } else {
          str = glue::glue("{props$type}{bounds} {name}")
        }
      } else if (props$type == "vector") {
        str = glue::glue("vector[{props$rows}]{bounds} {name}")
      } else if (props$type == "matrix") {
        str = glue::glue("matrix[{props$rows}, {props$cols}]{bounds} {name}")
      } else if (props$type == "assignment") {
        str = glue::glue("{name} = {props$expr}")
      } else if (props$type == "dist") {
        str = glue("{name} ~ {props$dist}({paste0(props$args, collapse = ', ')})")
      } else if (props$type == "loop") {
        return(private$make_loop(props$var, props$seq, props$expr))
      } else {
        stop("Unrecognized data type ", props$type)
      }
      
      if ("assignment" %in% names(props)) {
        str = glue("{str} = {props$assignment}")
      }
      
      paste0(str, ";")
    },
    
    make_loop = function(var, seq, expr) {
      expr = private$make_line(expr$name, expr$props)
      glue(
        "for ({var} in {seq}) {{
          {expr}
        }}"
      )
    }

  )
)

split_lines = function(string) {
  unlist(strsplit(string, "\n", fixed = TRUE))
}

add_space = function(string, n = 2) {
  space = paste0(rep(" ", n), collapse = "")
  paste0(space, string)
}
# var = "n"
# seq = "1:N"
# expr = "x + y"
# cat(paste(glue("for ({var} in {seq}) {{"), add_space(expr), "}", sep = "\n"))

