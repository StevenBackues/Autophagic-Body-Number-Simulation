# This file contains all of the functions required for type two clumping

require(pracma)
require(tidyverse)

tolerance <- 0.001

# this function returns a list of calculated points along the edge and active body ***** 3D READY *****
find_edge_points <- function(inc_rad, radius, rho, theta, phi, vac_rad, ...) {
  if (rho < vac_rad - radius) {
    return(NULL)
  }
  
  inc_rho <- vac_rad - inc_rad
  
  numer_edge <- (inc_rho ^ 2) + (rho ^ 2) - ((inc_rad + radius) ^ 2)
  denom_edge <- 2 * inc_rho * rho
  delta <- acos(numer_edge / denom_edge)
  
  # general rotation matrices
  r_z <- matrix(c(cos(theta + (pi / 2)), -sin(theta + (pi / 2)), 0, 0,
                  sin(theta + (pi / 2)),  cos(theta + (pi / 2)), 0, 0,
                                      0,                      0, 1, 0,
                                      0,                      0, 0, 1), 
                ncol = 4, byrow = TRUE)
  
  r_x <- matrix(c(1,                   0,                    0, 0,
                  0, cos(phi + (pi / 2)), -sin(phi + (pi / 2)), 0,
                  0, sin(phi + (pi / 2)),  cos(phi + (pi / 2)), 0,
                  0,                    0,                   0, 1),
                ncol = 4, byrow = TRUE)
  
  # disgonal specific rotation matrices
  r_z_d1 <- matrix(c(cos(pi / 4), -sin(pi / 4), 0, 0,
                     sin(pi / 4),  cos(pi / 4), 0, 0,
                               0,            0, 1, 0,
                               0,            0, 0, 1), 
                   ncol = 4, byrow = TRUE)
  r_z_d2 <- matrix(c(cos(3 * pi / 4), -sin(3 * pi / 4), 0, 0,
                     sin(3 * pi / 4),  cos(3 * pi / 4), 0, 0,
                                   0,                0, 1, 0,
                                   0,                0, 0, 1), 
                   ncol = 4, byrow = TRUE)
  
  # rectangular vector for the placed body, this shouldnt need to exist with better way to call function
  old_x <- rho * sin(phi) * cos(theta)
  old_y <- rho * sin(phi) * sin(theta)
  old_z <- rho * cos(phi)
  b_vec <- c(old_x, old_y, old_z)
  b_1 <- c(b_vec, 1)
  
  # pts 1 and 2 vary only by phi, easy!
  phi_new <- phi + delta
  pt_1 <- data_frame(radius = inc_rad,
                     rho = inc_rho,
                     theta = theta,
                     phi = phi_new,
                     x = inc_rho * sin(phi_new) * cos(theta),
                     y = inc_rho * sin(phi_new) * sin(theta),
                     z = inc_rho * cos(phi_new))
  
  phi_new <- phi - delta
  pt_2 <- data_frame(radius = inc_rad,
                     rho = inc_rho,
                     theta = theta,
                     phi = phi_new,
                     x = inc_rho * sin(phi_new) * cos(theta),
                     y = inc_rho * sin(phi_new) * sin(theta),
                     z = inc_rho * cos(phi_new))
  
  # pts 3 and 4 vary only by theta, less easy!
  b_1_prime <- b_1 %*% r_z
  b_1_prime2 <- b_1_prime %*% r_x

  theta_1_prime2 <- atan2(b_1_prime2[1, 2], b_1_prime2[1, 1])
  # pt 3
  theta_new_prime2 <- theta_1_prime2 + delta
  
  x_new_prime2 <- inc_rho * cos(theta_new_prime2)
  y_new_prime2 <- inc_rho * sin(theta_new_prime2)
  
  b_new_prime2 <- c(x_new_prime2, y_new_prime2, 0, 1)
  b_new_prime <- b_new_prime2 %*% inv(r_x)
  b_new <- b_new_prime %*% inv(r_z)

  theta_new <- atan2(b_new[1, 2], b_new[1, 1])
  phi_new <- acos(b_new[1, 3] / Norm(c(b_new[1, 1], b_new[1, 2], b_new[1, 3])))
  
  pt_3 <- data_frame(radius = inc_rad,
                     rho = inc_rho,
                     theta = theta_new,
                     phi = phi_new,
                     x = inc_rho * sin(phi_new) * cos(theta_new),
                     y = inc_rho * sin(phi_new) * sin(theta_new),
                     z = inc_rho * cos(phi_new))
 
  # pt 4
  theta_new_prime2 <- theta_1_prime2 - delta
  
  x_new_prime2 <- inc_rho * cos(theta_new_prime2)
  y_new_prime2 <- inc_rho * sin(theta_new_prime2)
  
  b_new_prime2 <- c(x_new_prime2, y_new_prime2, 0, 1)
  b_new_prime <- b_new_prime2 %*% inv(r_x)
  b_new <- b_new_prime %*% inv(r_z)
  
  theta_new <- atan2(b_new[1, 2], b_new[1, 1])
  phi_new <- acos(b_new[1, 3] / Norm(c(b_new[1, 1], b_new[1, 2], b_new[1, 3])))
  
  pt_4 <- data_frame(radius = inc_rad,
                     rho = inc_rho,
                     theta = theta_new,
                     phi = phi_new,
                     x = inc_rho * sin(phi_new) * cos(theta_new),
                     y = inc_rho * sin(phi_new) * sin(theta_new),
                     z = inc_rho * cos(phi_new))
  
  # pts 5, 6 are diagonal, not that easy!
  b_1_prime3 <- b_1_prime2 %*% r_z_d1
  theta_1_prime3 <- atan2(b_1_prime3[1, 2], b_1_prime3[1, 1])
  # pt 5
  theta_new_prime3 <- theta_1_prime3 + delta
  
  x_new_prime3 <- inc_rho * cos(theta_new_prime3)
  y_new_prime3 <- inc_rho * sin(theta_new_prime3)
  
  b_new_prime3 <- c(x_new_prime3, y_new_prime3, 0, 1)
  b_new_prime2 <- b_new_prime3 %*% inv(r_z_d1)
  b_new_prime <- b_new_prime2 %*% inv(r_x)
  b_new <- b_new_prime %*% inv(r_z)
  
  theta_new <- atan2(b_new[1, 2], b_new[1, 1])
  phi_new <- acos(b_new[1, 3] / Norm(c(b_new[1, 1], b_new[1, 2], b_new[1, 3])))
  
  pt_5 <- data_frame(radius = inc_rad,
                     rho = inc_rho,
                     theta = theta_new,
                     phi = phi_new,
                     x = inc_rho * sin(phi_new) * cos(theta_new),
                     y = inc_rho * sin(phi_new) * sin(theta_new),
                     z = inc_rho * cos(phi_new))

  # pt 6
  theta_new_prime3 <- theta_1_prime3 - delta
  
  x_new_prime3 <- inc_rho * cos(theta_new_prime3)
  y_new_prime3 <- inc_rho * sin(theta_new_prime3)
  
  b_new_prime3 <- c(x_new_prime3, y_new_prime3, 0, 1)
  b_new_prime2 <- b_new_prime3 %*% inv(r_z_d1)
  b_new_prime <- b_new_prime2 %*% inv(r_x)
  b_new <- b_new_prime %*% inv(r_z)
  
  theta_new <- atan2(b_new[1, 2], b_new[1, 1])
  phi_new <- acos(b_new[1, 3] / Norm(c(b_new[1, 1], b_new[1, 2], b_new[1, 3])))
  
  pt_6 <- data_frame(radius = inc_rad,
                     rho = inc_rho,
                     theta = theta_new,
                     phi = phi_new,
                     x = inc_rho * sin(phi_new) * cos(theta_new),
                     y = inc_rho * sin(phi_new) * sin(theta_new),
                     z = inc_rho * cos(phi_new))

  # pts 7, 8 work use the last rotation but work on this plane's phi instead of theta
  b_1_prime3 <- b_1_prime2 %*% r_z_d2
  theta_1_prime3 <- atan2(b_1_prime3[1, 2], b_1_prime3[1, 1])
  # pt 7
  theta_new_prime3 <- theta_1_prime3 + delta
  
  x_new_prime3 <- inc_rho * cos(theta_new_prime3)
  y_new_prime3 <- inc_rho * sin(theta_new_prime3)
  
  b_new_prime3 <- c(x_new_prime3, y_new_prime3, 0, 1)
  b_new_prime2 <- b_new_prime3 %*% inv(r_z_d2)
  b_new_prime <- b_new_prime2 %*% inv(r_x)
  b_new <- b_new_prime %*% inv(r_z)
  
  theta_new <- atan2(b_new[1, 2], b_new[1, 1])
  phi_new <- acos(b_new[1, 3] / Norm(c(b_new[1, 1], b_new[1, 2], b_new[1, 3])))
  
  pt_7 <- data_frame(radius = inc_rad,
                     rho = inc_rho,
                     theta = theta_new,
                     phi = phi_new,
                     x = inc_rho * sin(phi_new) * cos(theta_new),
                     y = inc_rho * sin(phi_new) * sin(theta_new),
                     z = inc_rho * cos(phi_new))

  # pt 8
  theta_new_prime3 <- theta_1_prime3 - delta
  
  x_new_prime3 <- inc_rho * cos(theta_new_prime3)
  y_new_prime3 <- inc_rho * sin(theta_new_prime3)
  
  b_new_prime3 <- c(x_new_prime3, y_new_prime3, 0, 1)
  b_new_prime2 <- b_new_prime3 %*% inv(r_z_d2)
  b_new_prime <- b_new_prime2 %*% inv(r_x)
  b_new <- b_new_prime %*% inv(r_z)
  
  theta_new <- atan2(b_new[1, 2], b_new[1, 1])
  phi_new <- acos(b_new[1, 3] / Norm(c(b_new[1, 1], b_new[1, 2], b_new[1, 3])))
  
  pt_8 <- data_frame(radius = inc_rad,
                     rho = inc_rho,
                     theta = theta_new,
                     phi = phi_new,
                     x = inc_rho * sin(phi_new) * cos(theta_new),
                     y = inc_rho * sin(phi_new) * sin(theta_new),
                     z = inc_rho * cos(phi_new))
  
  return(list(pt_1, pt_2, pt_3, pt_4, pt_5, pt_6, pt_7, pt_8))
}

find_touching_points <- function(body_a, body_b, inc_rad, vac_rad, ...) { # ***** IN 3d *****
  # determine plane for first set of rotations
  b_1 <- c(body_a$x, body_a$y, body_a$z)
  b_2 <- c(body_b$x, body_b$y, body_b$z)
  n <- pracma::cross(b_1, b_2)
  
  # inverse trig functions sometimes require correction
  theta_n <- atan(n[2] / n[1])
  phi_n <- acos(n[3] / Norm(n))
  
  # rotation matrices for first set
  r_z <- matrix(c(cos(theta_n + (pi / 2)), -sin(theta_n + (pi / 2)), 0, 0,
                  sin(theta_n + (pi / 2)),  cos(theta_n + (pi / 2)), 0, 0,
                  0,                       0,                        1, 0,
                  0,                       0,                        0, 1), 
                ncol = 4, byrow = TRUE)
  
  # this test seems to make the first set of rotations work correctly
  if (n[1] >= tolerance){
    r_x <- matrix(c(1,          0,           0, 0,
                    0, cos(phi_n), -sin(phi_n), 0,
                    0, sin(phi_n),  cos(phi_n), 0,
                    0,          0,           0, 1),
                    ncol = 4, byrow = TRUE)
  } else {
    r_x <- matrix(c(1,          0,           0, 0,
                    0, cos(pi - phi_n), -sin(pi - phi_n), 0,
                    0, sin(pi - phi_n),  cos(pi - phi_n), 0,
                    0,          0,           0, 1),
                    ncol = 4, byrow = TRUE)
  }
  
  # rotation of body_a to determine the direction for delta
  b_1 <- c(body_a$x, body_a$y, body_a$z, 1)
  b_1_prime <- b_1 %*% r_z
  b_1_prime2 <- b_1_prime %*% r_x
  
  # rotation of body_b for the same reason
  b_2 <- c(body_b$x, body_b$y, body_b$z, 1)
  b_2_prime <- b_2 %*% r_z
  b_2_prime2 <- b_2_prime %*% r_x
  
  # find new theta value for angle on this plane
  theta_1_prime2 <- atan2(b_1_prime2[1, 2], b_1_prime2[1, 1])
  theta_2_prime2 <- atan2(b_2_prime2[1, 2], b_2_prime2[1, 1])
  
  # this makes sure we are adding delta in the correct "direction" -- toward the secondary body
  if (theta_1_prime2 > theta_2_prime2) {
    # switch the entire bodies so our second rotations work
    temp <- body_b
    body_b <- body_a
    body_a <- temp
    # switch rectangular bodies
    temp <- b_2
    b_2 <- b_1
    b_1 <- temp
    # switch the rotated body vectors
    temp <- b_2_prime2
    b_2_prime2 <- b_1_prime2
    b_1_prime2 <- temp
    # and thetas
    temp <- theta_2_prime2
    theta_2_prime2 <- theta_1_prime2
    theta_1_prime2 <- temp
  }
  
  # finding angle alpha
  numer_alph <- (inc_rad + body_a$radius) ^ 2 + (body_a$radius + body_b$radius) ^ 2 - (inc_rad + body_b$radius) ^ 2
  denom_alph <- 2 * (inc_rad + body_a$radius) * (body_a$radius + body_b$radius)
  alpha <- acos(numer_alph / denom_alph)
  # finding angle beta
  numer_beta <- body_a$rho ^ 2 + (body_a$radius + body_b$radius) ^ 2 - body_b$rho ^ 2
  denom_beta <- 2 * body_a$rho * (body_a$radius + body_b$radius)
  beta <- acos(numer_beta / denom_beta)
  # finding angle gamma
  .gamma <- beta - alpha
  # finding new_rho
  new_rho <- sqrt( (inc_rad + body_a$radius) ^ 2 + body_a$rho ^ 2 - 2 * (inc_rad + body_a$radius) * body_a$rho * cos(.gamma) )
  # finding delta, the actual important angle for later calculations
  numer_touch <- body_a$rho ^ 2 + new_rho ^ 2 - (body_a$radius + inc_rad) ^ 2
  denom_touch <- 2 * body_a$rho * new_rho
  delta <- acos(numer_touch / denom_touch)
  
  # check whether delta should be added or subtracted
  if ((alpha > beta & abs(theta_1_prime2 - theta_2_prime2) < pi)
    | (alpha < beta & abs(theta_1_prime2 - theta_2_prime2) > pi)) {
    theta_new_prime2 <- theta_1_prime2 - delta
  } else {
    theta_new_prime2 <- theta_1_prime2 + delta
  }
  
  # find coordinates for new body on this plane
  x_new_prime2 <- new_rho * cos(theta_new_prime2)
  y_new_prime2 <- new_rho * sin(theta_new_prime2)
  # rotate back to original coordinates
  new_prime2 <- c(x_new_prime2, y_new_prime2, 0, 1)
  new_prime <- new_prime2 %*% inv(r_x)
  new_rect <- new_prime %*% inv(r_z)
  
  # define second normal vector
  n2 <- c(body_b$x - body_a$x, body_b$y - body_a$y, body_b$z - body_a$z)
  phi_n2 <- acos(n2[3] / Norm(n2))
  theta_n2 <- atan2(n2[2], n2[1])
  
  # normalize it then scale it to the length to find touching point, p
  n_hat <- n2 / Norm(n2)
  p <- b_1[1:3] + (n_hat * body_a$radius)
  
  # find vector, h, that connects P and the new body, project it onto second plane, find circle radius
  h <- p - new_rect[1, 1:3]
  proj_h <- h - (dot(n2, h) / (Norm(n2) ^ 2) * n2)
  viable_radius <- Norm(proj_h)
  
  # find circle origin trig style
  n_h <- viable_radius / tan(alpha)
  circle_origin <- c(body_a$x, body_a$y, body_a$z) + (n_h * n_hat)
  
  # transformation matrices for second set
  t_circle <- matrix(c(1, 0, 0, 0,
                       0, 1, 0, 0,
                       0, 0, 1, 0,
                       -circle_origin[1], -circle_origin[2], -circle_origin[3], 1),
                     ncol = 4, byrow = TRUE)
  
  r_z2 <- matrix(c(cos(theta_n2 + (pi / 2)), -sin(theta_n2 + (pi / 2)), 0, 0,
                   sin(theta_n2 + (pi / 2)),  cos(theta_n2 + (pi / 2)), 0, 0,
                   0,                         0,                        1, 0,
                   0,                         0,                        0, 1), 
                ncol = 4, byrow = TRUE)
  
  r_x2 <- matrix(c(1,           0,            0, 0,
                   0, cos(phi_n2), -sin(phi_n2), 0,
                   0, sin(phi_n2),  cos(phi_n2), 0,
                   0,           0,            0, 1), 
                ncol = 4, byrow = TRUE)
  
  # viable point 1
  coord <- c(viable_radius, 0, 0, 1)
  pt_1 <- make_touching_point(coord, inc_rad, r_x2, r_z2, t_circle)
  
  # viable point 2
  coord <- c((sqrt(2) / 2) * viable_radius, (sqrt(2) / 2) * viable_radius, 0, 1)
  pt_2 <- make_touching_point(coord, inc_rad, r_x2, r_z2, t_circle)
  
  # viable point 3
  coord <- c(0, viable_radius, 0, 1)
  pt_3 <- make_touching_point(coord, inc_rad, r_x2, r_z2, t_circle)
  
  # viable point 4
  coord <- c((-sqrt(2) / 2) * viable_radius, (sqrt(2) / 2) * viable_radius, 0, 1)
  pt_4 <- make_touching_point(coord, inc_rad, r_x2, r_z2, t_circle)
  
  # viable point 5
  coord <- c(-viable_radius, 0, 0, 1)
  pt_5 <- make_touching_point(coord, inc_rad, r_x2, r_z2, t_circle)
  
  # viable point 6
  coord <- c((-sqrt(2) / 2) * viable_radius, (-sqrt(2) / 2) * viable_radius, 0, 1)
  pt_6 <- make_touching_point(coord, inc_rad, r_x2, r_z2, t_circle)
  
  # viable point 7
  coord <- c(0, -viable_radius, 0, 1)
  pt_7 <- make_touching_point(coord, inc_rad, r_x2, r_z2, t_circle)
  
  # viable point 8
  coord <- c((sqrt(2) / 2) * viable_radius, (-sqrt(2) / 2) * viable_radius, 0, 1)
  pt_8 <- make_touching_point(coord, inc_rad, r_x2, r_z2, t_circle)
  
  # if (is_touching(pt_1, body_b) & is_touching(pt_2, body_b) &
  #     is_touching(pt_3, body_b) & is_touching(pt_4, body_b) &
  #     is_touching(pt_5, body_b) & is_touching(pt_6, body_b) &
  #     is_touching(pt_7, body_b) & is_touching(pt_8, body_b) &
  #     is_touching(pt_1, body_a) & is_touching(pt_2, body_a) &
  #     is_touching(pt_3, body_a) & is_touching(pt_4, body_a) &
  #     is_touching(pt_5, body_a) & is_touching(pt_6, body_a) &
  #     is_touching(pt_7, body_a) & is_touching(pt_8, body_a)) {
  #   print("Touching Points: TRUE")
  # } else {
  #   print("Touching Points: FALSE")
  # }

  return(list(pt_1, pt_2, pt_3, pt_4, pt_5, pt_6, pt_7, pt_8))
}

# this function saves on typing in the find touching points function
make_touching_point <- function(vec, body_rad, r_1, r_2, t) {
  vp_3 <- vec
  vp_2 <- vp_3 %*% inv(r_1)
  vp_1 <- vp_2 %*% inv(r_2)
  vp <- vp_1 %*% inv(t)
  
  x <- vp[1, 1]
  y <- vp[1, 2]
  z <- vp[1, 3]
  
  rho <- sqrt(x^2 + y^2 + z^2)
  phi <- acos(z / rho)
  theta <- atan2(y, x)
  
  ret <- data_frame(radius = body_rad,
                    rho = rho,
                    theta = theta,
                    phi = phi,
                    x = x,
                    y = y,
                    z = z)
  return(ret)
}

# this function returns the distance between two bodies ***** ADDED Z COORD *****
distance_between <- function(body_a, body_b) {
  check_matrix <- matrix(c(body_a$x, body_a$y, body_a$z,
                           body_b$x, body_b$y, body_b$z),
                         nrow = 2,
                         byrow = TRUE)
  as.numeric(dist(check_matrix))
}

# this function returns a logical based on if two bodies are touching
is_touching <- function(body_a, body_b) {
  if (distance_between(body_a, body_b) >= body_a$radius + body_b$radius - tolerance
      & distance_between(body_a, body_b) < body_a$radius + body_b$radius + tolerance) {
    TRUE
  } else {
    FALSE
  }
}

# this function returns a logical based on whether or not two bodies are overlapping
is_overlapping <- function(body_a, body_b) {
  if (distance_between(body_a, body_b) < body_a$radius + body_b$radius - tolerance) {
    TRUE
  } else {
    FALSE
  }
}

# ***** MAIN *****
type_two_clustering <- function(inc_radii, vac_radius) {
  bodies <- data_frame(radius = double(), rho = double(), theta = double(), phi = double(),
                       x = double(), y = double(), z = double())
  
  # First body sent to vacuole edge randomly ***** ADDED PHI, Z *****
  radius <- inc_radii[1]
  rho <- vac_radius - radius
  theta <- runif(1) * 2 * pi
  phi <- runif(1) * pi
  x <- rho * sin(phi) * cos(theta)
  y <- rho * sin(phi) * sin(theta)
  z <- rho * cos(phi)
  add_list <- list(radius = radius, rho = rho, theta = theta, phi = phi, x = x, y = y, z = z)
  bodies <- bind_rows(bodies, add_list)
  
  print("1 body placed")
  
  for (inc in seq_along(inc_radii[-1])) {
    valid_points <- data_frame()
    
    print("Finding viable points...")
    
    #find points where new body touches vacuole edge and one body
    edge_points_list <- bodies %>% pmap(find_edge_points, inc_rad = inc_radii[inc], vac_rad = vac_radius) %>% 
      unlist(recursive = FALSE)
    
    # create df of all possible combinations of placed body indices
    body_idx_vector <- 1:nrow(bodies)
    body_idx_df <- expand.grid(body_idx_vector, body_idx_vector)
    all_combn <- body_idx_df[which(body_idx_df$Var1 > body_idx_df$Var2), ] %>%
      pmap(function(Var1, Var2) { list(bodies[Var1,], bodies[Var2,]) })
    
    # find the touching bodies
    keep_pairs <- which(map(all_combn, ~is_touching(.[[1]], .[[2]])) == TRUE)
    touching_pairs <- all_combn[keep_pairs]
    
    # calculate and return points where new body touches two old bodies
    touching_points_list <- touching_pairs %>%
      map(~find_touching_points(inc_rad = inc_radii[inc], vac_rad = vac_radius, body_a = .[[1]], body_b = .[[2]])) %>%
      unlist(recursive = FALSE)
    
    # status report
    print(paste("Viable Edge Points:", length(edge_points_list), " Viable Touching points:", length(touching_points_list)))
    viable_points <- data.frame(bind_rows(edge_points_list, touching_points_list))
    print(paste("Total Viable points:", nrow(viable_points)))
    print("Deleting invalid points...")
    
    # We must delete any point which extends too far
    viable_points <- viable_points[which(viable_points$rho + viable_points$radius < vac_radius + tolerance), ]
    
    # create list containing all possible combinations of viable points and bodies
    viable_idx_vector <- 1:nrow(viable_points)
    viable_vs_bodies <- expand.grid(viable_idx_vector, body_idx_vector) %>%
      pmap(function(Var1, Var2) { list(viable_points[Var1, ], bodies[Var2, ]) })
    
    # determine which points overlap an existing body, return their index
    overlap_idx <- which(map(viable_vs_bodies, ~is_overlapping(.[[1]], .[[2]])) == TRUE)
    
    # if points need to be deleted for overlaps
    if (length(overlap_idx > 0)) {
      valid_points <- viable_vs_bodies[overlap_idx] %>% map(1) %>% bind_rows() %>% setdiff(x = viable_points)
    } else {
      valid_points <- viable_points %>% drop_na() %>% unique()
    }
    
    print(paste("Valid Points:", nrow(valid_points)))
    
    if (nrow(valid_points) > 0) {
      bodies <- bind_rows(bodies, sample_n(valid_points, 1))
    } else {
      message("No valid points found. Terminating body placement.")
      return(NULL)
    }
    
    print(paste(inc + 1, "bodies placed"))
    
  }
  return(bodies)
}