#' ci_sim
#'
#' This function performs bootstrap simulation for finding the confidence intervals of extrapolated feasibility counts with weighted proportions.
#'
#' @param n_cohorts the number of cohorts, or buckets
#' @param n_den a vector containing denominator counts for each cohort
#' @param n_rev a vector containing the number of patients reviewed in each cohort
#' @param n_case a vector containing the number of hits in each cohort for the feasibility question of interest
#' @param N_sim the number of simulations
#'
#' @return a list containing a tibble with simulation output, a tibble with count percentiles, a tibble with summarized results intended for excel use, a summary tibble, and a histogram
#'
#' @importFrom magrittr "%>%"
#' @import dplyr
#' @import tibble
#'
#' @export
#'
#' @examples
ci_sim = function(n_cohorts, n_den, n_rev, n_case, N_sim) {

  # Find the size of cohort with the most reviewed cases
  max_rev = max(n_rev)

  # create matrix of cohort samples
  Y_sample = matrix(NA, max_rev, n_cohorts)
  for (i in 1:n_cohorts) {
    Y_sample[,i] = rep(c(0,1,NA), c(n_rev[i]-n_case[i], n_case[i], max_rev-n_rev[i]))
  }

  ### Bootstrap simulation ###

  # Create bootstrap loop
  Y_boot = matrix(NA, N_sim, n_cohorts)

  for(c in 1:n_cohorts) {

    B_case = rep(NA, N_sim)

    for(b in 1:N_sim) {

      Y.b <- Y_sample[sample(1:n_rev[c], replace=TRUE),c]
      B_case[b] <- sum(Y.b)
    }

    Y_boot[,c] = B_case
  }

  ### Weighted Calculations ###

  # calculate the weighted sums of cases in each cohort, for every repetition
  weighted_sums = matrix(NA, N_sim, n_cohorts)

  for(c in 1:n_cohorts) {

    for(b in 1:N_sim) {

      weighted_sums[b,c] = ((Y_boot[b,c] / n_rev[c]) * (n_den[c] - n_rev[c])) + n_case[c]
    }
  }

  # Add all weighted sums, find Frequency and Percentage for each trial
  total_sums = as.data.frame(weighted_sums) %>%
    mutate(Frequency = rowSums(.[1:n_cohorts]),
           Percentage = Frequency / sum(n_den)) %>%
    select(Frequency, Percentage)

  # combine bootstrap data
  sim_data = as_tibble(bind_cols(as.data.frame(Y_boot), total_sums))


  # Make percentile table
  percentiles = c(0.01, 0.025, 0.05, .10, .20, .25, .30, .40, .50,
                  .60, .70, .75, .80, .90, .95, .975, .99)

  per_table = as.data.frame(quantile(sim_data$Frequency, percentiles)) %>%
    rownames_to_column() %>%
    select(Percentile = rowname,
           Mean_Freq = 'quantile(sim_data$Frequency, percentiles)') %>%
    as_tibble()


  # histogram of frequencies
  hist = sim_data %>%
    ggplot2::ggplot(ggplot2::aes(x=Frequency)) + ggplot2::geom_histogram()


  # create result table
  results = tibble(
    mean = round(quantile(sim_data$Frequency, 0.5),0),
    mean_p = round(quantile(sim_data$Percentage, 0.5),3),
    ci_95 = paste("(",round(quantile(sim_data$Frequency, .025),0),", ", round(quantile(sim_data$Frequency, .975),0), ")", sep=""),
    ci_95_p = paste("(", round(quantile(sim_data$Percentage, .025),3), ", ", round(quantile(sim_data$Percentage, .975),3), ")",sep=""),
    prob_90 = round(quantile(sim_data$Frequency, .10),0),
    prob_90_p = round(quantile(sim_data$Percentage, .10),3),
    prob_75 = round(quantile(sim_data$Frequency, .25),0),
    prob_75_p = round(quantile(sim_data$Percentage, .25),3)
  )


  # create summary table
  mean_freq = paste(round(quantile(sim_data$Frequency, 0.5),0), " (",round((quantile(sim_data$Percentage, 0.5))*100,1),"%)",sep = "")
  ci_low = paste(round(quantile(sim_data$Frequency, .025),0), " (",round((quantile(sim_data$Percentage, .025))*100,1),"%)",sep = "")
  ci_high = paste(round(quantile(sim_data$Frequency, .975),0), " (",round((quantile(sim_data$Percentage, .975))*100,1),"%)",sep = "")
  ci_95 = paste(ci_low, "-", ci_high)
  prob_90 = paste(round(quantile(sim_data$Frequency, .10),0), " (",round((quantile(sim_data$Percentage, .10))*100,1),"%)",sep = "")
  prob_75 = paste(round(quantile(sim_data$Frequency, .25),0), " (",round((quantile(sim_data$Percentage, .25))*100,1),"%)",sep = "")

  freq_table = cbind(mean_freq, ci_95, prob_90, prob_75) %>%
    as_tibble() %>%
    mutate(ID = "Frequency") %>%
    select(ID,
           "Mean" = mean_freq,
           "95% CI" = ci_95,
           "90% probability N≥" = prob_90,
           "75% probability N≥" = prob_75)


  summary_table = as_tibble(reshape2::dcast(reshape2::melt(freq_table, id.vars = "ID"), variable ~ ID) %>%
                              select("Statistic" = variable, Frequency))



  # Output results
  output = list(results_table = results,
                simulation_output = sim_data,
                percentile_table = per_table,
                histogram = hist,
                summary_table = summary_table)

  return(output)

}
