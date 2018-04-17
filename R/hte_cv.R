#' @import dplyr
#' @import purrr
#' @import tidyr
#' @import magrittr
#' @import caret
#' @import Matching

# from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3079915/
c_stat = function(est_ranking, time, event, treatment, weight=1) {
	event_data = data.frame(est_ranking=est_ranking, time=time, event=event, treatment=treatment, weight=weight, dummy=0) %>%
		filter(event)
	ordered_data = inner_join(event_data, event_data, by="dummy") %>%
		filter(time.x < time.y)
	denominator = (ordered_data$weight.x)^2 %>% sum # the ^2 is essentially because we compare two individuals (see ref 19 in paper above)
	numerator = ordered_data %>%
		filter(est_ranking.x < est_ranking.y) %>%
		mutate(weight=weight.x^2) %>%
		pull(weight) %>%
		sum()
	return(numerator/denominator) # the higher (closer to 1) the better
}

surv_mse = function(est_survival, time, event, treatment, weight=1) {
	data.frame(est_survival=est_survival, time=time, weight=weight) %>%
		mutate(mse = weight*(est_survival-time)^2) %>%
		pull(mse) %>%
		mean(na.rm=T)
}

value = function(est_treatment, time, event, treatment, weight=1) {
	data.frame(est_treatment=est_treatment, treatment=treatment, time=time, event=event, weight=weight) %>%
		filter(est_treatment == treatment) %>% # the "lucky" individuals
		arrange(-time) %>%
		mutate(weighted_n_alive = cumsum(weight)) %>%
		arrange(time) %>%
		mutate(km_surv = cumprod(1-(event*weight)/weighted_n_alive)) %>%
		mutate(dt = lead(time) - time) %>%
		mutate(rest_mean_survival = km_surv*dt) %>%
		pull(rest_mean_survival) %>% 
		sum(na.rm=T) # an estimate of the restricted mean survival time of the lucky patients
}

true_value = function(est_treatment, treated_mean, control_mean) {
	conditional_mean = est_treatment*treated_mean + (1-est_treatment)*control_mean
	mean(conditional_mean)
}

random_metric = function(){ return(0)}
