################################################################################
# Purpose: To analyse P and Q events
#
################################################################################
# This script processes P and Q data to define events
################################################################################
# Authors: Oscar Baez-Villanueva and Mauricio Zambrano-Bigiarini 
################################################################################
# Started: 17-September-2024
# Updates: 
################################################################################

#' Time series formatting (P and Q)
#'
#' @param q_obs zoo time series of observed streamflow for a defined period
#' @param p_obs zoo time series of observed streamflow for a defined period
#' @param par_beta numeric representing the filter parameter. Default value 
#' is 0.925 as recommended by Arnold and Allen (1999)
#'
#' @return zoo object with streamflow ('Q_tot', 'Q_bf', 'Q_dr', 'Q_dr_bool') 
#' and precipitation data ('P_obs', 'P_bool')
#' @export
#'
#' @examples
sepparate_streamflow = function(q_obs, p_obs, par_beta = 0.925){
  
  # Sepparating the base flow with a digital filter
  baseflow = hydroTSM::baseflow(q_obs, na.fill = 'linear', beta = par_beta)
  
  # Calculating direct runnof
  q_direct     = q_obs - baseflow
  q_direct_bol = round(q_direct, 0)
  
  # Calculating boolean event identification for streamflow
  q_boolean = q_direct_bol
  q_boolean[q_boolean <  1] = 0
  q_boolean[q_boolean >= 1] = 1
  
  
  # Calculating boolean event identification for precipitation (thers. 1 mm)
  p_boolean = p_obs
  p_boolean[p_boolean >= 1] = 1
  p_boolean[p_boolean < 1] = 0
  
  # Creating zoo with results
  res = cbind(q_obs, baseflow, q_direct, q_boolean, p_obs, p_boolean)
  names(res) = c('Q_tot', 'Q_bf', 'Q_dr', 'Q_dr_bool', 'P_obs', 'P_bool')
  
  return(res)
}

#' Event identification
#'
#' @param boolean boolean time series that identify an event (0 = no event, 
#' 1 = event)
#' @param threshold values equal to, or lower than this threshold will not be 
#' counted as an event
#'
#' @return
#' @export
#'
#' @examples
identify_events = function(boolean, threshold = 1) {
  # Initialise the result vector with zeroes
  result <- rep(0, length(boolean))
  
  # Counter for the number of events
  event_counter = 0
  
  # Start and end indexes for consecutive 1s
  start = NULL
  
  # Loop through the boolean vector
  for (i in seq_along(boolean)) {
    if (boolean[i] == 1) {
      if (is.null(start)) {
        start <- i  # Mark the start of the event
      }
    } else {
      if (!is.null(start)) {
        # Check the length of the event
        if ((i - start) > threshold) {
          event_counter = event_counter + 1
          # Assign the event number to the result
          result[start:(i-1)] <- event_counter  
        }
        # Reset the start
        start <- NULL 
      }
    }
  }
  
  # Check for an event that continues until the end of the vector
  if (!is.null(start)) {
    if ((length(boolean) - start + 1) > threshold) {
      event_counter <- event_counter + 1
      result[start:length(boolean)] = event_counter
    }
  }
  
  return(result)
}


#' Identify the date of the center of mass of the event and its total volume
#'
#' @param observations zoo time series with the observations
#' @param events vector identifying the number of event for the 'observations' 
#' time series
#'
#' @return
#' @export
#'
#' @examples
centroid_events = function(observations, events){
  
  # Getting the total number of identified events
  no_events = max(events)
  
  # Iterating trough each event to calculate the center of mass of each event
  #   and the respective date
  cm_date = vol_tot = start_event = end_event = NA
  for(i in 1:no_events){
    
    # Identifying the event
    pos_event = which(events == i)
    
    # If there is a one-day event, the date of the center of mass will be that
    #   specific day
    if(length(pos_event) > 1){
      
      # Extracting the values for the specific event
      values = observations[pos_event]
      dates  =  zoo::index(observations)[pos_event]
      
      # Calculating center of mass of the event according to the rounded date
      numeric_dates  = as.numeric(dates - min(dates)) + 1
      pos_cm         = round(sum(numeric_dates * values) / sum(values), 0)
      start_event[i] = as.character(dates[1])
      end_event[i]   = as.character(dates[length(dates)])
      cm_date[i]     = as.character(dates[pos_cm])
      vol_tot[i]     = sum(values)
      
    } else {
      start_event[i] = end_event[i] = cm_date[i] = as.character(zoo::index(observations)[pos_event])
      vol_tot[i]     = as.numeric(observations[pos_event])
    }
    
  }
  
  # Building the results data frame with the No_event, total_volume, and cm date
  res = data.frame(No_event     = 1:no_events, 
                   Total_volume = vol_tot, 
                   CM_date      = cm_date,
                   Start_event  = start_event,
                   End_event    = end_event)
  
  return(res)
  
}

