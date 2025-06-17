# Add cost complexity measure to mod$frame
# ChatGPT 6/3/25
# 
label_as_integer <- function(v) {
  v[v != "<leaf>"] <- as.character(seq_len(sum(v != "<leaf>")))
 # v[v == "<leaf>"] <- "\U1F343"
  v
}

###

vector_colors <- function(x,
                          color_low = "#deebf7",  # light blue
                          color_high = "#08306b", # dark blue
                          default_color = "gray90") {
  stopifnot(is.numeric(x))
  
  is_valid <- is.finite(x)
  valid_x <- x[is_valid]
  
  # Initialize all as default
  colors <- rep(default_color, length(x))
  
  # Return default if no variation or valid input
  if (length(valid_x) == 0 || diff(range(valid_x)) == 0) {
    return(colors)
  }
  
  # Normalize to [0, 1]
  x_norm <- (valid_x - min(valid_x)) / diff(range(valid_x))
  
  # Create gradient function
  ramp <- grDevices::colorRamp(c(color_low, color_high), space = "Lab")
  rgb_mat <- ramp(x_norm)
  hex_colors <- grDevices::rgb(rgb_mat[,1], rgb_mat[,2], rgb_mat[,3], maxColorValue = 255)
  
  colors[is_valid] <- hex_colors
  return(colors)
}

vector_colors <- function(x,
                          color_low = "#FF896A",  
                          color_high = "#deebf7", 
                          default_color = "gray90") {
  stopifnot(is.numeric(x))
  
  is_valid <- is.finite(x)
  valid_x <- x[is_valid]
  
  # Initialize all as default
  colors <- rep(default_color, length(x))
  
  # Return default if no variation or valid input
  if (length(valid_x) == 0 || diff(range(valid_x)) == 0) {
    return(colors)
  }
  
  # Normalize to [0, 1]
  x_norm <- (valid_x - min(valid_x)) / diff(range(valid_x))
  
  # Create gradient function
  ramp <- grDevices::colorRamp(c(color_low, color_high), space = "Lab")
  rgb_mat <- ramp(x_norm)
  hex_colors <- grDevices::rgb(rgb_mat[,1], rgb_mat[,2], rgb_mat[,3], maxColorValue = 255)
  
  colors[is_valid] <- hex_colors
  return(colors)
}

#### 

categorical_colors <- function(x,
                               palette = NULL,
                               default_color = "gray90",
                               na_color = "gray90") {
  # Convert input to factor
  x <- as.factor(round(x, 3))
  levels_x <- levels(x)
  n_levels <- length(levels_x)
  
  # Choose or generate palette
  if (is.null(palette)) {
    # Use RColorBrewer if available, otherwise rainbow
    if (requireNamespace("RColorBrewer", quietly = TRUE) &&
        n_levels <= 12) {
      palette <- RColorBrewer::brewer.pal(n_levels, "Set3")
    } else {
      palette <- grDevices::rainbow(n_levels)
    }
  }
  
  # Safety check
  if (length(palette) < n_levels) {
    stop("Not enough colors in palette for number of levels")
  }
  
  # Assign colors
  colors <- rep(na_color, length(x))
  is_valid <- !is.na(x)
  colors[is_valid] <- palette[as.integer(x[is_valid])]
  
  return(colors)
}

# from rpart.plot
get.class.stats <- function(x)
{
  # columns of yval2 for e.g. a two-level response are: fitted n1 n2 prob1 prob2
  yval2 <- x$frame$yval2
  if(NCOL(yval2) < 5)
    stop0("is.class.response(x) yet frame$yval2 is not ",
          "a matrix with five or more columns")
  fitted <- yval2[, 1] # fitted level as an integer
  if(NCOL(yval2) %% 2 == 0) { # new style yval2?
    stopifnot(colnames(yval2)[length(colnames(yval2))] == "nodeprob")
    nlev <- (ncol(yval2) - 2) / 2
  } else # old style yval2
    nlev <- (ncol(yval2) - 1) / 2
  stopifnot(nlev > 1)
  stopifnot(floor(nlev) == nlev)
  n.per.lev <- yval2[, 1 + (1:nlev), drop=FALSE]
  # aug 2012: commented out the following check because it
  # incorrectly fails when cases weights are used in the rpart model
  # stopifnot(sum(n.per.lev[1,]) == ntotal) # sanity check
  prob.per.lev <- yval2[, 1 + nlev + (1:nlev), drop=FALSE]
  # dec 2012: loosened following check to allow for numerical error
  stopifnot(abs(sum(prob.per.lev[1,]) - 1) < 1e-8) # sanity check
  list(yval2=yval2,
       fitted=fitted,
       nlev=nlev,
       # Aug 2019: ylevels is necessary for the multiclass model legend when
       # the last level in the response is unused in the training data and thus
       # does't appear in yval2 e.g. see "unusedlev" in the rpart.plot tests
       ylevels=attr(x, "ylevels"),
       n.per.lev=n.per.lev,
       prob.per.lev=prob.per.lev)
}

annotate_gt <- function(mod) {
  frame <- mod$frame
  is_leaf <- frame$var == "<leaf>"
  
  dev <- frame$dev
  node_nums <- rownames(frame)  # use character node IDs
  node_index <- setNames(seq_along(node_nums), node_nums)
  
  # Correctly find children using rpart's node numbering
  get_kids <- function(i) {
    this_node <- as.numeric(node_nums[i])
    left <- as.character(this_node * 2)
    right <- as.character(this_node * 2 + 1)
    kids <- c()
    if (left %in% node_nums) kids <- c(kids, node_index[left])
    if (right %in% node_nums) kids <- c(kids, node_index[right])
    return(kids)
  }
  
  get_descendant_leaves <- function(i) {
    out <- c()
    explore <- c(i)
    while (length(explore) > 0) {
      current <- explore[1]
      explore <- explore[-1]
      kids <- get_kids(current)
      if (length(kids) == 0) {
        out <- c(out, current)
      } else {
        explore <- c(explore, kids)
      }
    }
    return(out[frame$var[out] == "<leaf>"])
  }
  n <- nrow(frame)
  gt <- rep(NA, n) 
  subdev <- rep(NA, n) 
  numsplits <- rep(NA, n) 
  for (i in seq_len(nrow(frame))) {
    if (!is_leaf[i]) {
      leaves <- get_descendant_leaves(i)
      if (length(leaves) >= 2) {
        gt[i] <- (dev[i] - sum(dev[leaves])) / (length(leaves) - 1)
        subdev[i] <- sum(dev[leaves])
        numsplits[i] <- length(leaves)
      } else {
        gt[i] <- Inf
      }
    }
  }
  frame$gt <- gt
  frame$subdev <- subdev
  frame$numsplits <- numsplits
  return(frame)
}

