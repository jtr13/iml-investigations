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


# labeling functions
# 
# labels each node with class labels, class counts
# and deviance (equal to smaller count)
# default colors 
# 
node.fun.class.dev <- function(x, labs, digits, varlen) {
  stats <- get.class.stats(x)
  paste(stats$ylevels[1], stats$ylevels[2], "\n", stats$n.per.lev[,1], "      ", stats$n.per.lev[,2], 
        "\nmisclassified: ", x$frame$dev)    
}


node.fun.gt.calc <- function(x, labs, digits, varlen) x$frame$label

node.label <- function(x, labs, digits, varlen) x$frame$label

split.order.dev <- function(x, labs, digits, varlen) {
  l <- x$frame$var
  l[l != "<leaf>"] <- as.character(seq_len(sum(l != "<leaf>")))
  l <- ifelse(l == "<leaf>", "leaf", paste0("split #", l))
  paste(l, "\ndev: ", x$frame$dev)
}


# create colors for use in rpart.plot
# order must match order of mod$frame
# 
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

# annotate with g(t)  (misclassification reduction / # of splits)
# 
# Correctly find children using rpart's node numbering
get_kids <- function(i, node_index) {
  this_node <- as.numeric(names(node_index)[i])
  left <- as.character(this_node * 2)
  right <- as.character(this_node * 2 + 1)
  kids <- c()
  if (left %in% names(node_index)) kids <- c(kids, node_index[left])
  if (right %in% names(node_index)) kids <- c(kids, node_index[right])
  return(kids)
}

get_descendant_leaves <- function(frame, i, node_index) {
  out <- c()
  explore <- c(i)
  while (length(explore) > 0) {
    current <- explore[1]
    explore <- explore[-1]
    kids <- get_kids(current, node_index)
    if (length(kids) == 0) {
      out <- c(out, current)
    } else {
      explore <- c(explore, kids)
    }
  }
  return(out[frame$var[out] == "<leaf>"])
}

 
annotate_gt <- function(mod) {
  frame <- mod$frame
  is_leaf <- frame$var == "<leaf>"
  
  dev <- frame$dev
  node_index <- setNames(seq_along(rownames(frame)), rownames(frame))
  
 n <- nrow(frame)
  gt <- rep(NA, n) 
  subdev <- rep(NA, n) 
  numsplits <- rep(NA, n) 
  for (i in seq_len(nrow(frame))) {
    if (!is_leaf[i]) {
      leaves <- get_descendant_leaves(frame, i, node_index)
      if (length(leaves) >= 2) {
        gt[i] <- (dev[i] - sum(dev[leaves])) / (length(leaves) - 1)
        subdev[i] <- sum(dev[leaves])
        numsplits[i] <- length(leaves)
      }
    }
  }
  frame$gt <- gt
  frame$subdev <- subdev
  frame$numsplits <- numsplits
  return(frame)
}

get_descendants <- function(node, all_nodes) {
  # Recursive function to collect all descendants
  children <- c(2 * node, 2 * node + 1)
  descendants <- children[children %in% all_nodes]
  for (child in descendants) {
    descendants <- c(descendants, get_descendants(child, all_nodes))
  }
  return(descendants)
}


gt_calc_colors <- function(mod, showpruned = TRUE) {
  if (nrow(mod$frame) == 1) {
    "#00A86B50"
  } else {
    mingt <- min(mod$frame$gt, na.rm = TRUE)
    mingt_nodes <- na.omit(mod$frame$nn[mod$frame$gt == mingt])
    mingt_node_cp <- mod$frame$complexity[mod$frame$nn == mingt_nodes][1]
    to_be_pruned <- unique(unlist(
      lapply(mingt_nodes, get_descendants, all_nodes = mod$frame$nn)
    ))
    # c6e1ff is used in classification_pruning.qmd, be careful!
    if (showpruned) {
      dplyr::case_when(
        mod$frame$nn %in% to_be_pruned ~ "grey90",
        mod$frame$gt == mingt ~ "#c6c9ff",
        !is.na(mod$frame$gt) ~ "#c6e1ff",
        .default = "#00A86B50"
      )
      } else {
      dplyr::case_when(
        !is.na(mod$frame$gt) ~ "#c6e1ff",
        .default = "#00A86B50" 
      )
      }
  }
}

create_label <- function(mod, gt = FALSE, cp = TRUE,
                         showmin = TRUE) {

    # deal with 1 node situation
  if (nrow(mod$frame) == 1) {
    return(paste("   leaf   \ndev: ", mod$frame$dev))
  }
  
  # create leaf labels
  label <- ifelse(mod$frame$var == "<leaf>",
    paste("   leaf   \ndev: ", mod$frame$dev), "")
  
  if (gt) {
    label_gt <- ifelse(
      mod$frame$var != "<leaf>",
      paste0("g(t) = (", mod$frame$dev, "-",
             mod$frame$subdev, ") / \n (",
             mod$frame$numsplits, "-1) = ",
             round(mod$frame$gt, 2)), "")
      label <- paste0(label, label_gt, "\n")
    }
  
    if (cp) {
      label_cp <- ifelse(mod$frame$var != "<leaf>",
        paste0("cp: ",
        round(mod$frame$gt/mod$frame$dev[1], 5)), "")
      label <- paste0(label, label_cp, "\n")
    }
  
    if (showmin) {
    mingt <- min(mod$frame$gt, na.rm = TRUE)  
    label_showmin <- dplyr::case_when(
      mod$frame$gt == mingt ~ "*minimum",
      .default = "")
    label <- paste0(label, label_showmin)
    }
  sub("\\n+$", "", label)
  }

