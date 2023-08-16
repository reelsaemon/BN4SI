# Preliminaries -----------------------------------------------------------

gc()

# extended functionality
require(lubridate)
require(data.table)
require(tidyr)
require(plyr)
require(dplyr)
require(stringr)
require(rlist)
require(DescTools)

# model packages
require(bnlearn)

# evaluation metrics
require(caret)
require(Metrics)
require(StatRank)
require(comparator)

# parallelization
require(foreach)
require(doParallel)

source('src/BN4SI_functions.R')


# Model setup -------------------------------------------------------------

set.seed(26477)

# constants
totalCores <- detectCores()
usedCores <- totalCores-1
n_samples <- 10000
training_size <- 2/3
test_size <- 1-training_size
horizon <- 1

# Data import -------------------------------------------------------------

event_log_data <- as.data.table(read.csv("data/sim_event_log.csv", header = T, sep = ","))
id_colname <- 'order_id'
activity_colname <- 'station'
timestamp_colname <- 'timestamp_in'
event_log_data <- event_log_data[, c(id_colname, activity_colname, timestamp_colname), with=FALSE]

# Data preprocessing ------------------------------------------------------

event_log_data[, (activity_colname) := as.character(get(activity_colname))]
event_log_data[, trace_length := .N, by=get(id_colname)]
max_trace_length <- event_log_data[, max(trace_length)]
variants <- event_log_data[, .(variants=paste0(get(activity_colname), collapse="-")), by=get(id_colname)]
event_log_data[, FirstTimestamp := min(get(timestamp_colname)), by=get(id_colname)]
event_log_data <- event_log_data[order(FirstTimestamp), ]

# chronological train test split
trace_ids <- event_log_data[, unique(get(id_colname))]
train_ids <- trace_ids[1:floor((training_size)*length(trace_ids))]
test_ids <- trace_ids[(floor((training_size)*length(trace_ids))+1):length(trace_ids)]

# # random train test split
# train_ids <- sample(event_log_data[, unique(get(id_colname))], floor(training_size*event_log_data[, length(unique(get(id_colname)))]), replace=FALSE)
# test_ids <- setdiff(event_log_data[, unique(get(id_colname))], train_ids)

training_data <- event_log_data[get(id_colname) %in% train_ids, ]
max_trace_length_train <- training_data[, max(trace_length)]
activities_train <- training_data[, unique(get(activity_colname))]
test_data <- event_log_data[get(id_colname) %in% test_ids, ]

training_data_variants <- training_data[, .(variants=paste0(get(activity_colname), collapse="-")), by=get(id_colname)]
test_data_variants <- test_data[, .(variants=paste0(get(activity_colname), collapse="-")), by=get(id_colname)]


# place-annotation and preprocessing of datasets
training_data_wide <- make_wide_format(training_data, activity_colname, id_colname, binary=TRUE, preserve_exec=FALSE)
test_data_wide <- make_wide_format(test_data, activity_colname, id_colname, binary=TRUE, preserve_exec=TRUE, input_max_trace_length=max_trace_length_train, input_activities=activities_train)
test_data_wide_all_info <- copy(test_data_wide)
test_data_wide_all_info_cols <- setdiff(colnames(test_data_wide_all_info), c(id_colname, "execution_order"))
test_data_wide <- test_data_wide[, c(id_colname, "execution_order", setdiff(colnames(training_data_wide), id_colname)), with=FALSE] # limit information in test data set


# Network construction and fitting ----------------------------------------

network_structure <- make_sparse_binary_structure_check_max_parents(training_data, activity_colname, id_colname, horizon=horizon)
fitted_network <- bn.fit(network_structure, training_data_wide[, nodes(network_structure), with=FALSE])

# Next Activity/Remaining Trace Prediction --------------------------------

times <- paste0("_at_time_", str_pad(string=1:(max_trace_length+1), width=1+floor(log10(max_trace_length+1)), pad="0", side="left"))

close(file("parallel_log.txt"), open="w")
cluster <- makeCluster(usedCores, outfile="parallel_log.txt")
registerDoParallel(cluster)
pb <- txtProgressBar(min=0, max=nrow(test_data_wide), initial=0, width=50, style=3)

predicted_traces_list <- foreach(i = 1:nrow(test_data_wide), .packages=c("data.table", "bnlearn", "tidyr", "stringr")) %dopar% {

  evidence <- as.list(test_data_wide[i, lapply(.SD, function(x) if(any(is.na(x))){NULL} else {x}), .SDcols=nodes(network_structure)])
  evidence_nodes <- names(evidence)
  
  queried_times <- paste0("_at_time_", str_pad(string=(test_data_wide[i, execution_order]+1):(max_trace_length+1), width=1+floor(log10(max_trace_length+1)), pad="0", side="left"))
  query_nodes <- setdiff(nodes(network_structure), evidence_nodes)
  
  evidence_times <- paste0("_at_time_", str_pad(string=(test_data_wide[i, execution_order]-(horizon-1)):(test_data_wide[i, execution_order]), width=1+floor(log10(max_trace_length+1)), pad="0", side="left"))
  evidence_nodes <- evidence_nodes[grep(paste0(evidence_times, collapse="|"), evidence_nodes)]
  evidence <- evidence[evidence_nodes]
  
  if(length(query_nodes) > 0) {
    
    query_nodes_list <- lapply(queried_times, function(x) query_nodes[grep(x, query_nodes)])
    names(query_nodes_list) <- queried_times
    
    predicted_trace <- character(0)
    list_index <- 1
    current_prediction <- "placeholder"
    
    while(!grepl("^act_finished", current_prediction) & (list_index >= 1 & list_index <= length(query_nodes_list))) {
      
      current_query_nodes <- query_nodes_list[[list_index]]
      
      current_node_query <- suppressWarnings(clean_multi_cpdist(fit=fitted_network, current_query_nodes, evidence=evidence, n_obs=10000))
      
      if(is.null(current_node_query)) {
        
        print("we could not produce a query...")
        
        current_prediction <- NA
        # names(current_prediction) <- names(query_nodes_list)[list_index]
        predicted_trace <- c(predicted_trace, current_prediction)
        
        # predicted_trace <- as.vector(predicted_trace)
        predicted_activity <- predicted_trace[1]
        
        setTxtProgressBar(pb, i)
        return(list(predicted_trace=predicted_trace, predicted_activity=predicted_activity))
      }
      
      occ_probs <- unlist(lapply(current_node_query[names(current_node_query) %in% current_query_nodes], function(x) x[[TRUE]]))
      
      # handle unseen data
      if(!(any(occ_probs!=0))|any(is.nan(occ_probs))) {
        
        current_prediction <- NA
        predicted_trace <- c(predicted_trace, current_prediction)
        
        predicted_activity <- predicted_trace[1]
        
        setTxtProgressBar(pb, i)
        return(list(predicted_trace=predicted_trace, predicted_activity=predicted_activity))
        
      } else {
        
        current_prediction <- str_remove(names(occ_probs[which.max(occ_probs)]), names(query_nodes_list)[list_index])
        
        if(length(current_prediction) > 1) {
          
          print("we have multiple most probable activities!")
        }
        
        names(current_prediction) <- names(query_nodes_list)[list_index]
        predicted_trace <- c(predicted_trace, current_prediction)
        
      }
      
      list_index <- list_index + 1
    }
    
    
    if(predicted_trace[length(predicted_trace)] == "act_finished") {
      
      predicted_trace <- as.vector(predicted_trace[1:(length(predicted_trace)-1)])
      predicted_activity <- predicted_trace[1]
      
      if(predicted_trace=="act_finished") {
        
        predicted_trace <- character(0)
        predicted_activity <- character(0)
      }
    } else {
      
      print("this should never happen...")
      predicted_trace <- as.vector(predicted_trace)
    }
    
  } else {
    
    predicted_trace <- character(0)
    predicted_activity <- character(0)
  }
  
  setTxtProgressBar(pb, i)
  
  return(list(predicted_trace=predicted_trace, predicted_activity=predicted_activity))
}

close(pb)
stopCluster(cluster)

eval_frame <- copy(test_data_wide)

eval_frame[, predicted_trace := lapply(predicted_traces_list, function(x) if(!any(is.na(x))) {x[['predicted_trace']][x[['predicted_trace']]!="act_finished"]} else {NA})]

eval_frame[, predicted_activity := lapply(predicted_traces_list, function(x) if(!any(is.na(x))) {x[['predicted_activity']][x[['predicted_activity']]!="act_finished"]} else {NA})]
eval_frame[, predicted_activity := as.character(predicted_activity)]
eval_frame[, predicted_activity := ifelse(predicted_activity=="NA", NA, predicted_activity)]

known_activities <- apply(test_data_wide_all_info[, test_data_wide_all_info_cols, with=FALSE], 1, function(x) names(x[!is.na(x)][x[!is.na(x)]==TRUE]))
known_activities <- lapply(known_activities, function(x) list(str_remove(x, paste0(times, collapse="|"))))

test_data_wide_all_info[, (test_data_wide_all_info_cols) := lapply(.SD, function(x) unique(na.omit(x))), .SDcols=test_data_wide_all_info_cols, by=get(id_colname)]

all_activities <- apply(test_data_wide_all_info[, test_data_wide_all_info_cols, with=FALSE], 1, function(x) list(names(x[x==TRUE])))
all_activities <- lapply(all_activities, unlist)
all_activities <- lapply(all_activities, function(x) list(str_remove(x, paste0(times, collapse="|"))))

remaining_activities <- lapply(1:length(known_activities), function(x) if(test_data_wide[x, execution_order] < length(all_activities[[x]][[1]])) {all_activities[[x]][[1]][(length(known_activities[[x]][[1]])+1):length(all_activities[[x]][[1]])]} else {character(0)})

eval_frame[, complete_trace := lapply(all_activities, function(x) str_remove_all(paste0(unlist(x), collapse="-"), "-act_finished"))]

eval_frame[, actual_trace := lapply(remaining_activities, function(x) x[x!="act_finished"])]

eval_frame[, actual_activity := lapply(remaining_activities, function(x) x[1][x[1]!="act_finished"])]
eval_frame[, actual_activity := as.character(actual_activity)]

# evaluation with accuracy and similarity metrics

eval_frame[, current_activity := lag(actual_activity, 1), by=get(id_colname)]
eval_frame <- eval_frame[current_activity!="act_finished" | is.na(current_activity), ]

eval_frame[, last_index:=last(execution_order), by=get(id_colname)]
eval_frame <- eval_frame[execution_order!=last_index, ]

# next activity prediction accuracy
print(paste0("Next Activity Accuracy: ", eval_frame[!is.na(predicted_activity), sum(predicted_activity == actual_activity)]/nrow(eval_frame)))

# remaining trace similarity
similarities <- c()
for(index in 1:nrow(eval_frame)) {
  
  similarities <- c(similarities, manual_ndls(eval_frame[index, actual_trace], eval_frame[index, predicted_trace]))
}
similarities[is.na(similarities)] <- 0

print(paste0("Suffix NDLS: ", mean(similarities))) # this is our similarity
