# generates vacuoles  with log-normal distributed radii with inputed vmu and vsigma 
# randomly positions them on an axis, and slices them with "slice 1," including only those giving a minimum observed radius (vacmin)
# For each sliced vacuole, slice 2 then uses function genBalls to fill it with bodies, where body number is normally distributed (bm and bsd) and body size is lognormal (mu and sigma), and CompactSpheres clumps them 
# Slice2 outputs the max radius captured in the slice for each ball in the slice
# then discards everything below the recognition limit (reclim)
# repeats this for process for the input number of cells, combining all the results
# outputs the input data as well as the observed radii of the body crossections and the number of crossections per vacuole 
# repeats all of this for the input number of times (repeats), varying mean body number, and outputs all of the results to a text file "mean size and number.txt"
resim = function (mu = 5.07, sigma = 0.339, vmu = 6.95, vsigma = 0.25, bm = 30, 
                  bsd = 4, reclim = 45, vacmin = 300, cells = 200, repeats = 4,
                  iter_parameter = "bm", clustering = 2, cores = 3) {
  
  require(MultiRNG)
  require(tidyverse)
  require(foreach)
  require(doParallel)
  source("Body Number Analysis Plugin.R")
  source("type_two_clumping.R")
  source("glpk_clumping.R")
  
  # sets timestamp for writing to files
  start_time <- Sys.time()
  char_date <- as.Date(start_time)
  
  if (file.exists("input_functions.R")) {
    source("input_functions.R")
    # accepts user input, assigns inputs to correct variables
    input <- read_user_input(mu, sigma, vmu, vsigma, bm, bsd, reclim, vacmin, cells, repeats, iter_parameter, clustering, cores)
    
    mu <- input[[1]]
    sigma <- input[[2]]
    vmu <- input[[3]]
    vsigma <- input[[4]]
    bm <- input[[5]]
    bsd <- input[[6]]
    reclim <- input[[7]]
    vacmin <- input[[8]]
    cells <- input[[9]]
    repeats <- input[[10]]
    iter_parameter <- input[[11]]
    clustering <- input[[12]]
    cores <- input[[13]]
  }
  
  # set up cluster for parallel processing
  cl <- makeCluster(cores)
  registerDoParallel(cl)
  
  input_record <- data_frame(parameter = c("mu", "sigma", "vmu", "vsigma", "bm",
                                           "bsd", "reclim", "vacmin", "cells",
                                           "repeats", "iter_parameter", "clustering",
                                           "cores"),
                             value = c(mu, sigma, vmu, vsigma, bm, bsd, reclim,
                                       vacmin, cells, repeats, iter_parameter,
                                       clustering, cores))
  
  # needed for parallel iteration
  start_bm <- bm
  start_bsd <- bsd
  
  # creates a lognormal distribution of vacuoles, based on vmu and vsigma, then filters it to make sure that none are smaller than the vacuole recognition limit (vacmin)
  RRALL = NULL
  VACS = NULL
  vacuole_data = data_frame()
  body_summary = data_frame(input_mean_rad = double(), input_sd_rad = double(), obs_mean_rad = double(),
                     obs_sd_rad = double(), input_mean_num = integer(), input_sd_num = integer(), obs_per_vac = double())
  
  vacuole_data <- foreach (i = 1:repeats, 
                           .combine = "rbind",
                           .export = ls(envir=globalenv()),
                           .packages = c("tidyverse", "MultiRNG", "pracma")) %dopar% {
    # make sure iterations work correctly
    if (iter_parameter == "bsd") {
      bsd <- start_bsd + i
      bm <- start_bm
    } else if (iter_parameter == "bm") {
      bm <- start_bm + i
      bsd <- start_bsd
    } else {
      print(paste("Invalid iterate parameter:", iter_parameter))
      print("Check inputs")
      break
    }
    
    # makes vacuoles, following a lognormal size distribution
    vacradl = rlnorm(10*cells, mean=vmu, sd=vsigma)
    # filters the vacuoles to discard any smaller than the vacuole recognition limit
    vacradf = vacradl[vacradl > vacmin]
    write(vacradf, file = "original vacuoles.txt", sep = "\t", ncolumns=1)
    # Thickness of a slice, T (based on thickness of the ultrathin sections used for TEM)
    T = 70 
    # slices the vacuoles
    s1 = slice1(vacradf, vacmin, cells, T)
    s1$ivacuoles = s1$ivacuoles[, 1:cells]
    # creates bodies in each vacuole and calculates areas of those bodies caught in slice
    s2 = slice2(s1$ivacuoles, mu, sigma, bm, bsd, reclim, T, clustering)
    RRALL = c(RRALL, s2$RRall)
    meanrad = mean(RRALL)
    SDrad = sd(RRALL)
    num = length(RRALL)
    numc = (num/ncol(s1$ivacuoles))
    write(RRALL, file = "body radii.txt", sep = "\t", ncol=1)

    # build a data frame made up of some general summary data for the body distribution
    outpt = data_frame(input_mean_rad = mu,
                       input_sd_rad = sigma, 
                       obs_mean_rad = meanrad,
                       obs_sd_rad = SDrad,
                       input_body_mean = bm, 
                       input_body_sd = bsd,
                       obs_per_vac = numc)
    body_summary = bind_rows(body_summary, outpt)

    # build vacuole data as a dataframe
    stats = data_frame(input_body_mean = rep.int(bm, nrow(s2$vacuoles)), input_body_sd = rep.int(bsd, nrow(s2$vacuoles)))
    VACS = data_frame(true_vac_rad = s2$vacuoles[, 1], 
                      sliced_vac_rad = s2$vacuoles[, 2], 
                      true_bodies = s2$vacuoles[, 3],
                      sliced_bodies = s2$vacuoles[, 4])
    VACUOLES = bind_cols(stats, VACS)
    return(VACUOLES)
  }
  
  # end multicore processing
  stopCluster(cl)
  
  # output
  sorted_ks <- run_ks_tests(vacuole_data, char_date)
  plot_split_violins(vacuole_data, sorted_ks, char_date)
  plot_QQs(vacuole_data, char_date)
  plot_means(vacuole_data, char_date)
  write_data(vacuole_data, char_date)
  write_data(body_summary, char_date)
  write_data(input_record, char_date)
  
  print("done")
  print(Sys.time() - start_time)
}

# Takes a slice of thickness T of a box of vacuoles
slice1 = function(vacuoles, vacmin, cells, T){
  # Sets the size of the box to 10 times the mean vacuole radius 
  cellSize = mean(vacuoles)*10
  # randomly positions the vacuoles along the Z axis in this box
  vacpos = runif(length(vacuoles))* cellSize
  #takes a random slice
  z = 0.5*cellSize
  RVM = NULL
    for (j in 1:length(vacuoles)) {
     # the distance between the slice and the center of the vacuole
    zoffsetv = (vacpos [j] - z)
    # if the center of the vacuole is in the slice
    if (abs(zoffsetv) < T/2) {
      rv = vacuoles[j] 
      rvm = matrix (c(zoffsetv, rv, vacuoles[j]), nrow=3, ncol=1) 
      RVM = cbind (RVM, rvm)
      # if the vacuole is cut by the edge of the slice
    } else if ((abs(zoffsetv)-T/2) < vacuoles[j]) {
      # the crossectional radius of the vacuole at the edge of the slice
      rv = sqrt(vacuoles[j]^2-(abs(zoffsetv)-T/2)^2) 
      # verifies that the sliced radius is greater than the vacuole recognition limit
      if (rv >= vacmin) {
         rvm = matrix (c(zoffsetv, rv, vacuoles[j]), nrow=3, ncol=1) 
         RVM = cbind (RVM, rvm)
      } else NULL
    } else NULL
  }
  # ivacuoles is a matrix of vacuoles included in the slice, including the distance between their center and the slice, the sliced radius and the original radius.
  write.csv (RVM, file="ivacuoles.csv")
  list (ivacuoles = RVM)
}

# Slice 2
# Takes a slice (thickness, T) of each included vacuole, at the appropriate distance from the center of the vacuole
slice2 = function(ivacuoles, mu, sigma, bm, bsd, reclim, T, clustering){
  RRall = NULL
  SL2 = NULL
  for (v in 1:ncol(ivacuoles)) {
    rr = NULL
    RR = NULL
    print ("vacuole")
    print (v)
    # determines how many bodies for this vacuole - the number follows a normal distribution with input mean and sd
    bods = rnorm (1, mean=bm, sd=bsd)
    # rounds the number of bodies to a whole number
    bodies = round(bods)
    # if there is at least one body, will take a slice
    if (bodies >=1){
      # if there are at least two bodies, runs the clumping routine
      if (bodies >=2){
        # creates bodies with sizes following a lognormal distribution with input mu and sigma
        r=rlnorm(bodies,mean=mu,sd=sigma)
        vacrad = ivacuoles [3, v]
        cellvol = (2*vacrad)^3
        # for each vacuole, generates bodies with a lognormal distribution based on mu and sigma, and positions and clusters them using genBalls
        # balls is now a matrix of the position (in three dimensions) and radius of each ball
        if (clustering == 0 | clustering == 1) {
          balls = tryCatch ( {
            ret = genBalls(r,vacrad, clustering)
            Rm = matrix (ret$R, nrow=1, ncol=bodies)
            rbind(Rm, ret$compact)
          }, error = function(e) {
            message("Something went wrong in body generation and positioning, so this vacuole has been ommitted from analysis")
            message(paste("Original error:", e))
            return(NULL)
          } )
        } else if (clustering == 2) {
          balls = tryCatch ( {
            ret = type_two_clustering(r,vacrad)
            Rm = matrix (ret$radius, nrow=1, byrow = TRUE)
            trans_x <- ret$x + vacrad
            trans_y <- ret$y + vacrad
            trans_z <- ret$z + vacrad
            rbind(Rm, matrix(c(trans_x, trans_y, trans_z), nrow=3, byrow = TRUE))
          }, error = function(e) {
            message("Something went wrong in body generation and positioning, so this vacuole has been ommitted from analysis")
            message(paste("Original error:", e))
            return(NULL)
          } )
        }
        # this is the case there is just one body, so no clustering is done, and "balls" is a matrix with the size and position of that one body
      } else {
        r=rlnorm(1,mean=mu,sd=sigma);
        vacrad = ivacuoles [3, v]
        X = generate.point.in.sphere (1,3)
        scale=vacrad
        # generate random distances from the center, with maximum = vacrad
        D = runif(1)*scale
        # scale each point by a random distance, so that all are inside vacuole (sphere), and recenter that sphere to positive coordinates
        pos = X*D+scale
        # transpose matrix, so that rows are coordinates and columns are points, as required later
        tpos = t(pos)
        Rm = matrix (r, nrow=1, ncol=1)
        balls = rbind(Rm, tpos)
      }
      # verifies that "balls" isn't null (as would be generated if the total volume of bodies was larger than that of the vacuole.  If it is, outputs "NAs" 
      if (is.null(balls)) {
        RRout = NA
        ballnum = NA
        nb = NA
        # now slices the bodies
      } else {
        # slice location pulled from slice 1 outputs, expressed relative to the bottom of the current vacuole
        z = vacrad - ivacuoles [1, v]
        # loops for the number of balls
        for (j in 1:ncol(balls)) {
          # the distance between the slice and the center of the ball
          zoffset = abs(balls[4,j] - z)
          # the center of the ball is in the slice, so observed radius = radius of the ball
          if (zoffset < T/2) {
            rr = balls[1, j] 
            # in this case, the ball is cut by the edge of the slice, so the observed radius is calculated
          } else if ((zoffset-T/2) < balls[1, j]) {
            rr = sqrt(balls[1, j]^2-(zoffset-T/2)^2) 
            # if the ball isn't in the slice
          } else rr=0
          RR = c(RR,rr)
        }
        RRout = RR[RR>reclim]
        # Counts the number of ball slices that were observed in this vacuole (nb)
        nb = length(RRout)
        # Counts the original number of balls that were generated for the vacuole (ballnum)
        ballnum = ncol(balls)
      }
      #if there are no bodies
    } else { RRout=NULL
    ballnum = 0
    nb = 0
    }
    # record of how many observed bodies per slice, with relevant information
    sl2 = matrix (c(ivacuoles[3, v], ivacuoles[2, v], ballnum, nb), ncol=4, nrow=1)
    SL2 = rbind (SL2, sl2)
    # record of the all of the observed body radii
    RRall = c(RRall, RRout)
  }
  list (RRall=RRall, vacuoles=SL2)
}
