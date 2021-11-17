predict_iteration <- function(object, new_data, iter) {

  check_serialization(object)
  # stopifnot(object$pred_type == "regression")
  stopifnot(class(new_data) == "data.frame")

  ## Checking stuff Missingness
	if (!object$use_missing_data){
		nrow_before = nrow(new_data)
		new_data = na.omit(new_data)
		if (nrow_before > nrow(new_data)){
			cat(nrow_before - nrow(new_data), "rows omitted due to missing data. Try using the missing data feature in \"build_bart_machine\" to be able to predict on all observations.\n")
		}
	}
	if (nrow(new_data) == 0){
		stop("No rows to predict.\n")
	}

  java_bart_machine <- object$java_bart_machine
  n                 <- nrow(new_data)
  new_data          <- pre_process_new_data(new_data, object)

  # if(object$pred_type == "regression") {
  y_hat <- .jcall(java_bart_machine, "[D", "getSingleGibbsSamplePrediction", 
                  .jarray(new_data, dispatch = TRUE), as.integer(iter - 1), 
                  .jevalArray)
  # }

  return(y_hat)

}
