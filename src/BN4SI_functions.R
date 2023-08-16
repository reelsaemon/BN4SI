# generic follows graph function
get_follows_graph <- function(log, activity_colname, id_colname, plot=TRUE, add_text=TRUE, sort=TRUE, horizon=1, only_plot=FALSE, transition_threshold=1) {
  
  # params
  # horizon     if horizon=1 a directly follows graph is computed
  #             if 1 < horizon <= n a follows graph with n steps into the future trace is computed. If horizon>=n with n being the maximum trace length in the event log
  #             we get a eventually follows graph for all possible traces and trace lengths
  # 
  # sort=TRUE   sorts the activities with respect to the most followers, i.e. activities that have many followers are assigned a high spot in the matrix
  #             with eventually follows graphs those activities are most likely the ones earlier in the process
  #             with directly follows graphs those activities are most likely choice node activities, i.e. activities after which many different activities can be carried out
  
  activities <- sort(unlist(unique(log[[activity_colname]])))
  n_activities <- length(activities)
  
  ids <- sort(unlist(unique(log[[id_colname]])))
  
  follows_matrix <- matrix(0, nrow=n_activities, ncol=n_activities, dimnames=list(activities, activities))
  
  for(current_id in ids) {
    
    # print(current_id)
    
    current_trace_indices <- which(log[[id_colname]] == current_id)
    current_trace <- log[current_trace_indices, ][[activity_colname]]
    trace_length <- length(current_trace)
    
    i <- 1
    
    while(i <= trace_length) {
      
      current_activity <- current_trace[i]
      current_remaining_trace <- current_trace[-(1:i)][1:min(horizon, length(current_trace) - i)]
      
      if(!any(is.na(current_remaining_trace))) {
        for(act in current_remaining_trace) {
          
          follows_matrix[as.character(current_activity), as.character(act)] <- follows_matrix[as.character(current_activity), as.character(act)] + 1
        }
      }
      
      i <- i+1
    }
  }
  
  # apply threshold
  follows_matrix[apply(follows_matrix, 2, function(x) x/log[, table(get(activity_colname))]) < (1-transition_threshold)] <- 0
  
  if(sort==TRUE) {
    
    non_zero_entries <- follows_matrix
    non_zero_entries[non_zero_entries > 0] <- 1
    non_zero_followers <- apply(non_zero_entries, 1, sum)
    
    sorted_order <- names(sort(non_zero_followers, decreasing=TRUE))
    
    # reorder matrix rows/columns
    follows_matrix <- follows_matrix[sorted_order, sorted_order]
  }
  
  if(plot==TRUE) {
    
    image(t(follows_matrix[nrow(follows_matrix):1,]), xlab="to", ylab="from", axes=FALSE, col=heat.colors(dim(follows_matrix)[1]**2), useRaster=TRUE, main=paste0("Follows Graph of ", deparse(substitute(log)), " with horizon ", horizon))
    abline(1, -1)
    axis(1, at=seq(0,1,1/(ncol(follows_matrix)-1)), labels=rownames(follows_matrix))
    axis(2, at=seq(0,1,1/(ncol(follows_matrix)-1)), labels=rev(colnames(follows_matrix)))
    
    # values as text
    if(add_text==TRUE) {
      
      text(expand.grid(seq(0,1,1/(ncol(follows_matrix)-1)), seq(0,1,1/(ncol(follows_matrix)-1))), labels=t(follows_matrix[nrow(follows_matrix):1,]), cex=0.66)
    }
    
  }
  
  if(plot & only_plot) {
    return(invisible(NULL))
  } else {
  
    return(follows_matrix)
  }
  
}

# clean query output from the Process-aware Bayesian Network
clean_multi_cpdist <- function(fit, nodes, evidence, method='lw', n_obs, consider_NA_samples=FALSE) {
  
  dist <- tryCatch(cpdist(fit=fit, nodes=nodes, evidence=evidence, method=method, n=n_obs), error=function(cond) {return(NULL)})
  
  clean_output <- as.list(rep(NA, length(nodes)))
  names(clean_output) <- nodes
  
  if(is.null(dist)) {
    
    return(NULL)
  }
  
  if(consider_NA_samples) {
    for(node in nodes) {
      clean_output[[node]] <- round(table(dist[, node], useNA='always')/n_obs, digits=4)
    }
    return(clean_output)
  } else {
    for(node in nodes) {
      clean_output[[node]] <- round(table(dist[, node])/length(na.omit(dist[, node])), digits=4)
    }
    return(clean_output)
  }
}


# function for footprint matrix from follows graph
footprint <- function(follows_graph) {
  
  fg_dimnames <- dimnames(follows_graph)
  
  # initialize matrix with same dimensions as a directly/eventually follows graph
  footprint_matrix <- matrix(NA, nrow=nrow(follows_graph), ncol=ncol(follows_graph), dimnames=fg_dimnames)
  
  # check the follows graph for existing connections/transitions and fill them into the footprint matrix
  for(row in 1:nrow(footprint_matrix)) {
    
    for(col in row:ncol(footprint_matrix)) {
      
      if(follows_graph[row, col] > 0 & follows_graph[col, row] == 0) {
        
        footprint_matrix[row, col] <- ">"
        footprint_matrix[col, row] <- "<"
      } else if(follows_graph[row, col] == 0 & follows_graph[col, row] > 0) {
        
        footprint_matrix[row, col] <- "<"
        footprint_matrix[col, row] <- ">"
      } else if(follows_graph[row, col] > 0 & follows_graph[col, row] > 0){
        
        footprint_matrix[row, col] <- "||"
        footprint_matrix[col, row] <- "||" 
      } else {
        
        footprint_matrix[row, col] <- "-"
        footprint_matrix[col, row] <- "-"
      }
    }
  }
  
  return(footprint_matrix)
}


# function for automatically converting event log to wide format
make_wide_format <- function(log_data, activity_colname, id_colname, binary=TRUE, preserve_exec=FALSE, input_max_trace_length=NULL, input_activities=NULL) {
  
  require(mltools)
  require(tidyr)
  
  log_data <- copy(log_data[, c(id_colname, activity_colname), with=FALSE])
  original_max_trace_length <- max(log_data[, .N, by=get(id_colname)][, N])
  
  if(!is.null(input_max_trace_length)) {
    
    original_max_trace_length <- input_max_trace_length
  }
  
  log_data <- log_data[, c(get(activity_colname), rep("finished", max(0, original_max_trace_length - .N + 1))), by=get(id_colname)]
  colnames(log_data) <- c(id_colname, activity_colname)
  log_data[, time_slice := paste0("at_time_", str_pad(string=1:.N, width=1+floor(log10(.N)), pad="0", side="left")), by=get(id_colname)]
  
  if(preserve_exec) {
    
    log_data[, execution_order := 1:.N, by=get(id_colname)]
    
    
    if(!is.null(input_activities)) {
      
      own_activities <- log_data[get(activity_colname)!="finished", sort(unique(get(activity_colname)))]
      log_data[, (activity_colname) := factor(get(activity_colname), levels=c(union(own_activities, input_activities), "finished"))]
    } else {
      
      log_data[, (activity_colname) := factor(get(activity_colname))]
    }
    
    # first widening step
    wide_log_data <- as.data.table(pivot_wider(log_data[, c(id_colname, activity_colname, "time_slice", "execution_order"), with=FALSE], names_from="time_slice", values_from=activity_colname, names_expand=TRUE))
    
    if(!is.null(input_activities)) {
      
      activities <- union(own_activities, input_activities)
    } else {
      
      activities <- log_data[get(activity_colname)!="finished", sort(unique(get(activity_colname)))]
    }
    
    
    slice_nodes <- log_data[, unique(time_slice)]
    
    # carry values for each column after they are known for the first time
    wide_log_data[, (slice_nodes) := lapply(.SD, carry_helper), .SDcols=slice_nodes, by=get(id_colname)]
  } else {
    
    # first widening step
    wide_log_data <- as.data.table(pivot_wider(log_data[, c(id_colname, activity_colname, "time_slice"), with=FALSE], names_from="time_slice", values_from=activity_colname))
    
    if(!is.null(input_activities)) {
      
      activities <- union(own_activities, input_activities)
    } else {
      
      activities <- log_data[get(activity_colname)!="finished", sort(unique(get(activity_colname)))]
    }
  }
  
  time_columns <- sort(paste0("at_time_", str_pad(string=1:(original_max_trace_length + 1), width=1+floor(log10(original_max_trace_length + 1)), pad="0", side="left")))
  
  activity_at_time_columns <- apply(expand.grid(paste0("act_", c(activities, "finished"), "_"), time_columns), 1, paste, collapse="")
  
  # factorize activity columns
  wide_log_data[, (time_columns) := lapply(.SD, function(x) factor(x, levels=c(activities, "finished"))), .SDcols=time_columns]
  
  if(binary) {
    
    # second widening step
    
    # perform one hot encoding and recoding to boolean values
    one_hot_slices <- mltools::one_hot(wide_log_data[, time_columns, with=FALSE])
    one_hot_slices[, (colnames(one_hot_slices)) := lapply(.SD, function(x) factor(as.logical(x), levels=c(TRUE, FALSE))), .SDcols=colnames(one_hot_slices)]
    
    # add one hot encoded boolean values to original data table
    wide_log_data[, (activity_at_time_columns) := one_hot_slices]
    
    # remove time columns
    wide_log_data[, (time_columns) := NULL]
  }
  
  return(wide_log_data)
}

# function for automatically converting event log to wide format for prefix prediction task
make_wide_format_prefix <- function(log_data, activity_colname, id_colname, binary=TRUE, preserve_exec=FALSE, input_max_trace_length=NULL, input_activities=NULL) {
  
  require(mltools)
  require(tidyr)
  
  log_data <- copy(log_data[, c(id_colname, activity_colname), with=FALSE])
  original_max_trace_length <- max(log_data[, .N, by=get(id_colname)][, N])
  
  if(!is.null(input_max_trace_length)) {
    
    original_max_trace_length <- input_max_trace_length
  }
  
  log_data <- log_data[, c(get(activity_colname), rep("finished", max(original_max_trace_length - .N + 1))), by=get(id_colname)]
  colnames(log_data) <- c(id_colname, activity_colname)
  log_data[, time_slice := paste0("at_time_", str_pad(string=1:.N, width=1+floor(log10(.N)), pad="0", side="left")), by=get(id_colname)]
  
  if(preserve_exec) {
    
    log_data[, execution_order := 1:.N, by=get(id_colname)]
    
    
    if(!is.null(input_activities)) {
      
      own_activities <- log_data[get(activity_colname)!="finished", sort(unique(get(activity_colname)))]
      log_data[, (activity_colname) := factor(get(activity_colname), levels=c(union(own_activities, input_activities), "finished"))]
    } else {
      
      log_data[, (activity_colname) := factor(get(activity_colname))]
    }
    
    # first widening step
    wide_log_data <- as.data.table(pivot_wider(log_data[, c(id_colname, activity_colname, "time_slice", "execution_order"), with=FALSE], names_from="time_slice", values_from=activity_colname, names_expand=TRUE))
    
    if(!is.null(input_activities)) {
      
      activities <- union(own_activities, input_activities)
    } else {
      
      activities <- log_data[get(activity_colname)!="finished", sort(unique(get(activity_colname)))]
    }
    
    
    slice_nodes <- log_data[, unique(time_slice)]
    
    # carry values for each column after they are known for the first time
    wide_log_data[, (slice_nodes) := lapply(.SD, carry_helper_prefix), .SDcols=slice_nodes, by=get(id_colname)]
  } else {
    
    # first widening step
    wide_log_data <- as.data.table(pivot_wider(log_data[, c(id_colname, activity_colname, "time_slice"), with=FALSE], names_from="time_slice", values_from=activity_colname))
    
    if(!is.null(input_activities)) {
      
      activities <- union(own_activities, input_activities)
    } else {
      
      activities <- log_data[get(activity_colname)!="finished", sort(unique(get(activity_colname)))]
    }
  }
  
  time_columns <- sort(paste0("at_time_", str_pad(string=1:(original_max_trace_length + 1), width=1+floor(log10(original_max_trace_length + 1)), pad="0", side="left")))
  
  activity_at_time_columns <- apply(expand.grid(paste0("act_", c(activities, "finished"), "_"), time_columns), 1, paste, collapse="")
  
  # factorize activity columns
  wide_log_data[, (time_columns) := lapply(.SD, function(x) factor(x, levels=c(activities, "finished"))), .SDcols=time_columns]
  
  if(binary) {
    
    # second widening step
    
    # perform one hot encoding and recoding to boolean values
    one_hot_slices <- mltools::one_hot(wide_log_data[, time_columns, with=FALSE])
    one_hot_slices[, (colnames(one_hot_slices)) := lapply(.SD, function(x) factor(as.logical(x), levels=c(TRUE, FALSE))), .SDcols=colnames(one_hot_slices)]
    
    # add one hot encoded boolean values to original data table
    wide_log_data[, (activity_at_time_columns) := one_hot_slices]
    
    # remove time columns
    wide_log_data[, (time_columns) := NULL]
  }
  
  return(wide_log_data)
}

# helper function to carry the first non-NA element of a vector onto the remaining elements of the vector
carry_helper <- function(vector, value_to_carry=NULL) {
  
  if(length(vector)==1) {
    
    return(vector)
  }
  
  
  current_element <- NA
  
  for(i in 1:(length(vector)-1)) {
    
    current_element <- vector[i]
    
    if(!is.na(current_element)) {
        
      if(is.null(value_to_carry)) {
        
        vector[i+1] <- current_element
      } else {
        
        vector[i+1] <- value_to_carry
      }
    }
  }
  
  return(vector)
}

# helper function to carry the first non-NA element of a vector onto the remaining elements of the vector
carry_helper_prefix <- function(vector) {
  
  if(length(vector)==1) {
    
    return(vector)
  }
  
  current_element <- NA
  
  for(i in length(vector):1) {
    current_element <- vector[i]
    if(!is.na(current_element)) {
      vector[i-1] <- current_element
    }
  }
  
  return(vector)
}

# function for automatic generation of a sparse binary approach BN structure with checking the max parent count
make_sparse_binary_structure_check_max_parents <- function(log, activity_colname, id_colname, horizon, keep_non_occurring=TRUE, max_parent_count=29) {
  
  log <- copy(log[, c(id_colname, activity_colname), with=FALSE])
  max_trace_length <- max(log[, .N, by=get(id_colname)][, N])
  
  # number of parameters in the binary matrix approach is 2^(n_parents+1)
  # take 29 because 2^31 is exactly INT_MAX
  max_parent_count <- max_parent_count
  
  # append finished activity at the end of traces
  log <- log[, c(get(activity_colname), rep("finished", max_trace_length - .N + 1)), by=get(id_colname)]
  colnames(log) <- c(id_colname, activity_colname)
  
  log[, place := 1:.N, by=get(id_colname)]
  
  log[, act_at_time := paste0("act_", get(activity_colname), "_at_time_", str_pad(string=place, width=1+floor(log10(max(place))), side="left", pad = "0"))]
  
  follows_graph <- get_follows_graph(log, activity_colname="act_at_time", id_colname=id_colname, plot=FALSE, add_text=FALSE, sort=FALSE, horizon=horizon)
  
  footprint <- footprint(follows_graph)
  
  occurring_transitions <- data.frame(from=character(0), to=character(0))
  for(row in dimnames(footprint)[[1]]) {
    
    for(col in dimnames(footprint)[[2]]) {
      
      if(footprint[row, col] == ">") {
        occurring_transitions <- rbind.data.frame(occurring_transitions,data.frame(from=row, to=col))
      }
    }
  }
  
  # to be able to still get non-NA predictions out of the network for activities at process slices that they never occured on in the training data
  # we still include nodes for all activities at all slices
  # only the edges are sparsely distributed across the graph
  if(keep_non_occurring) {
    
    activities <- log[, unique(get(activity_colname))]
    time_columns <- paste0("at_time_", str_pad(string=1:(max_trace_length+1), width=1+floor(log10(max_trace_length+1)), pad="0", side="left"))
    activity_at_time_nodes <- apply(expand.grid(paste0("act_", activities, "_"), time_columns), 1, paste, collapse="")
    sparse_network <- cpdag(empty.graph(nodes=activity_at_time_nodes))
    # adjust arcs of the network to the ones in occurring_transitions
    arcs(sparse_network) <- occurring_transitions
  } else {
    
    activity_at_time_nodes <- unique(c(occurring_transitions[, "from"], occurring_transitions[, "to"]))
    sparse_network <- cpdag(empty.graph(nodes=activity_at_time_nodes))
    # adjust arcs of the network to the ones in occurring_transitions
    arcs(sparse_network) <- occurring_transitions
  }
  
  # remove parents if parent count exceeds max_parent_count
  for(node in nodes(sparse_network)) {
    
    # check if parent count exceeds max_parent_count
    if(length(sparse_network$nodes[[node]]$parents) > max_parent_count) {
      
      print(node)
      
      parent_activities <- sparse_network$nodes[[node]]$parents
      
      parent_connections_to_node <- list()
      
      for(pa in parent_activities) {
        
        parent_connections_to_node[[pa]] <- follows_graph[pa, node]
      }
      
      most_influential_parents <- sort(unlist(parent_connections_to_node), decreasing=TRUE)[1:max_parent_count]
      
      sparse_network$nodes[[node]]$parents <- names(most_influential_parents)
    }
  }
  
  
  return(sparse_network)
}

manual_ndls <- function(actual_trace, predicted_trace) {
  
  # actual_trace and predicted_trace are the sequences as lists of character vectors
  
  if(is.null(predicted_trace[[1]])) {
    
    return(NA)
  } else {
    
    if(length(actual_trace[[1]])==0 & length(predicted_trace[[1]])==0) {
      
      ndls <- 1
      return(ndls)
    }
    
    distance <- DamerauLevenshtein()(actual_trace, predicted_trace)
    normalized_distance <- distance/max(length(actual_trace[[1]]), length(predicted_trace[[1]]))
    ndls <- 1-normalized_distance
    return(ndls)
  }
}