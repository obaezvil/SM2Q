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
  intensity = duration = NA
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
      intensity[i]   = max(values)
      duration[i]    = length(dates)
      
    } else {
      start_event[i] = end_event[i] = cm_date[i] = as.character(zoo::index(observations)[pos_event])
      vol_tot[i]     = as.numeric(observations[pos_event])
      intensity[i]   = as.numeric(observations[pos_event])
      duration[i]    = 1
    }
    
  }
  
  # Building the results data frame with the No_event, total_volume, and cm date
  res = data.frame(No_event     = 1:no_events, 
                   Total_volume = vol_tot, 
                   Intensity    = intensity,
                   Duration     = duration,
                   CM_date      = cm_date,
                   Start_event  = start_event,
                   End_event    = end_event)
  
  return(res)
  
}


#' Attribute Q events to defined P events 
#'
#' @param p_centroid_events data.frame resulting from the 'centroid_events' 
#' function applied to P events
#' @param q_centroid_events data.frame resulting from the 'centroid_events' 
#' function applied to Q events
#'
#' @return This function retrieves a data.frame with the respective lag time,
#' between P and Q events and the respective information of both events 
#' (intensity, duration; and starting, ending, and CM times)
#' @export
#'
#' @examples
match_P2Q_events = function(p_centroid_events, q_centroid_events){
  
  # Identifying the number of P events
  no_p_events = nrow(p_centroid_events)
  
  # Iterating trough number of P events
  days_apart = detected_event = NA
  for(i in 1:no_p_events){
    
    cm_p_event  =  as.Date(p_centroid_events$CM_date[i])
    cm_q_events = as.Date(q_centroid_events$CM_date)
    
    # Find forward dates (dates greater or equal than cm_p_event)
    forward_dates = cm_q_events[cm_q_events >= cm_p_event]
    
    if(length(forward_dates) == 0){
      
      days_apart[i]      = -999
      detected_event[i]  = -999
      
    } else {
      
      # Find the closest forward date and computing the distance in days
      closest_forward_date = min(forward_dates)
      days_apart[i]        = as.numeric(closest_forward_date - cm_p_event)
      detected_event[i]    = which(as.Date(q_centroid_events$CM_date) %in% closest_forward_date)
      
      # Checking distance from previous event (case that the p event contributes 
      #   to an existing event)
      if(i > 1){
        previous_event       = detected_event[i] - 1
        previous_cm          = q_centroid_events$CM_date[previous_event]
        days_apart_prevEvent = as.numeric(as.Date(cm_p_event) - as.Date(previous_cm))
        
        # If previous event is closer to the P centroid than the forward event and
        #   the CM of the P event lies within the previous Q event's start and 
        #   ending dates, the P event is attributed to the previous event
        if(days_apart_prevEvent < days_apart[i]){
          
          q_start_date = q_centroid_events$Start_event[previous_event]
          q_end_date   = q_centroid_events$End_event[previous_event]
          
          # Check wether the p event lies within the start and end dates of the
          #   previous event
          is_within_range <- as.Date(cm_p_event) >= as.Date(q_start_date) & 
            as.Date(cm_p_event) <= as.Date(q_end_date)
          if(is_within_range){
            
            days_apart[i]        = days_apart_prevEvent
            detected_event[i]    = previous_event
            
          }
          
        }
      } # END if (for the first event) i different to 1
      
    } # END first if
    
  } # END for
  
  # Constructing the resulting object
  res = data.frame(No_P_event      = 1:no_p_events,
                   Related_Q_event = detected_event,
                   Days_apart      = days_apart,
                   P_volume        = p_centroid_events$Total_volume,
                   P_intensity     = p_centroid_events$Intensity,
                   P_duration      = p_centroid_events$Duration,
                   P_start         = p_centroid_events$Start_event,
                   P_end           = p_centroid_events$End_event,
                   P_CM            = p_centroid_events$CM_date,
                   Q_volume        = NA,
                   Q_peak          = NA,
                   Q_duration      = NA,
                   Q_start         = NA,
                   Q_end           = NA,
                   Q_CM            = NA)
  
  # Iterative process to add Q event properties
  for(i in 1:no_p_events){
    
    # find corresponding Q event
    qevent_pos = which(q_centroid_events$No_event %in% res$Related_Q_event[i])
    qevent     = q_centroid_events[qevent_pos,]
    
    if(nrow(qevent) != 0){
      
      res$Q_volume[i]        = qevent$Total_volume
      res$Q_peak[i]          = qevent$Intensity
      res$Q_duration[i]      = qevent$Duration
      res$Q_start[i]         = qevent$Start_event
      res$Q_end[i]           = qevent$End_event
      res$Q_CM[i]            = qevent$CM_date
      
    } else {
      
      res$Q_volume[i] =  res$Q_peak[i] =  res$Q_duration[i] = -999
      res$Q_start[i]  =   res$Q_end[i] = res$Q_CM[i]        = -999
      
    }
    
  } # END for
  
  return(res)
  
}

