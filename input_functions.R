# This scipt contains the input functions for the resim function
# if not present, resim will take arguments trough the console when called
# or as defaults

# the two following functions make it cleaner to write our general console input reading
numeric_user_in <- function(parameter) {
  while (TRUE) {
    inp <- readline(prompt = paste("Enter value for ", parameter, ": ", sep = ""))
    if (grepl("[0-9]+\\.?[0-9]*", inp)) {
      return(as.numeric(inp))
    } else {
      warning("Value must be numeric", call. = FALSE, immediate. = TRUE)
      next
    }
  }
}

allowed_user_in <- function(parameter, allowed) {
  while (TRUE) {
    inp <- readline(prompt = paste("Enter value for ", parameter, ": ", sep = ""))
    if (inp %in% allowed) {
      return(as.character(inp))
    } else {
      warning(paste("Value must be in", allowed), call. = FALSE, immediate. = TRUE)
      next
    }
  }
}

# this function will allow users to easily select input parameters
read_user_input <- function(mu, sigma, vmu, vsigma, bm, bsd, reclim, vacmin,
                            cells, repeats, iter_parameter, clustering, cores) {
  inputs <- list(length = 13)
  
  while (TRUE) {
    
    # how should we get input values?
    inp_type <- as.numeric(readline(prompt = paste("Please choose (1, 2, or 3) from the following list of inputs:\n",
                                                   "\t1) Read input from function defaults\n",
                                                   "\t2) Read input from 'resim_input.txt'\n",
                                                   "\t3) Enter input directly through console\n", sep = "")))
    
    if (inp_type == 1) {
      # inputs passed directly through resim call
      inputs[[1]] <- mu
      inputs[[2]] <- sigma
      inputs[[3]] <- vmu
      inputs[[4]] <- vsigma
      inputs[[5]] <- bm
      inputs[[6]] <- bsd
      inputs[[7]] <- reclim
      inputs[[8]] <- vacmin
      inputs[[9]] <- cells
      inputs[[10]] <- repeats
      inputs[[11]] <- iter_parameter
      inputs[[12]] <- clustering
      inputs[[13]] <- cores
    } else if (inp_type == 2) {
      # read inputs from file
      inputs <- tryCatch({
        file_input <- read_lines("resim_input.txt")
        input_clean <- sapply("^.*: (.*)$", sub, "\\1", file_input)
        defaults <- list(length = 13)
        
        defaults[[1]] <- as.double(input_clean[1])
        defaults[[2]] <- as.double(input_clean[2])
        defaults[[3]] <- as.double(input_clean[3])
        defaults[[4]] <- as.double(input_clean[4])
        defaults[[5]] <- as.integer(input_clean[5])
        defaults[[6]] <- as.integer(input_clean[6])
        defaults[[7]] <- as.numeric(input_clean[7])
        defaults[[8]] <- as.numeric(input_clean[8])
        defaults[[9]] <- as.integer(input_clean[9])
        defaults[[10]] <- as.integer(input_clean[10])
        
        if (input_clean[11] %in% c("bm", "bsd")) {
          defaults[[11]] <- input_clean[11]
        } else {
          defaults[[11]] <- NA
        }
        
        if (input_clean[12] %in% c("0", "1", "2")) {
          defaults[[12]] <- as.integer(input_clean[12])
        } else {
          defaults[[12]] <- NA
        }
        
        if (as.integer(input_clean[13]) %in% c(1:detectCores())) {
          defaults[[13]] <- as.integer(input_clean[13])
        } else {
          defaults[[13]] <- NA
        }
        
        return(defaults)
      }, error = function(e) {
        warning("File missing or corrupt. Using function defaults", call. = FALSE, immediate. = TRUE)
        print(paste("Original", e))
        next
      })
    } else if (inp_type == 3) {
      # read directly from console
      inputs[[1]] <- numeric_user_in("mu of lognormal distribution of body radii")
      inputs[[2]] <- numeric_user_in("sigma of lognormal distribution of body radii")
      inputs[[3]] <- numeric_user_in("mu of lognormal distribution of vacuole radii")
      inputs[[4]] <- numeric_user_in("sigma of lognormal distribution of vacuole radii")
      inputs[[5]] <- numeric_user_in("mean of normal distribution of body number")
      inputs[[6]] <- numeric_user_in("standard deviation of normal distribution of body number")
      inputs[[7]] <- numeric_user_in("minimum observable body cross-section radius")
      inputs[[8]] <- numeric_user_in("minimum observable vacuole cross-section radius")
      inputs[[9]] <- numeric_user_in("simulated sample size")
      inputs[[10]] <- numeric_user_in("number of iterating runs")
      inputs[[11]] <- allowed_user_in("parameter to iterate over, bm for body mean, bsd for body std. dev.", c("bm", "bsd"))
      inputs[[12]] <- as.integer(allowed_user_in("clustering method, 0: none, 1: glpk, 2: type two", c("0", "1", "2")))
      print(paste("Your system has", detectCores(), "CPU cores"))
      inputs[[13]] <- as.integer(allowed_user_in("number of cores for parallel processing", c(1:detectCores())))
    } else {
      print("Please enter either 1, 2, or 3 to make your selection")
      next
    }
    
    return(inputs)
  }
}
