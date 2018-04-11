#' @import dplyr
#' @import purrr
#' @import tidyr
#' @import magrittr
#' @import caret
#' @import Matching

# from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3079915/
c_stat = function(est_rel_risk, time, event, treatment, IPCW=1) {
	event_data = data.frame(est_rel_risk=est_rel_risk, time=time, event=event, treatment=treatment, IPCW=IPCW, dummy=0) %>%
		filter(event)
	# print(event_data)
	ordered_data = inner_join(event_data, event_data, by="dummy") %>%
		filter(time.x < time.y)
		# print(ordered_data)
	denominator = (ordered_data$IPCW.x)^2 %>% sum
	numerator = ordered_data %>%
		filter(est_rel_risk.x > est_rel_risk.y) %>%
		mutate(IPCW=IPCW.x^2) %>%
		pull(IPCW) %>%
		sum()
	return(numerator/denominator) # the higher (closer to 1) the better
}

value = function(est_effect, treatment, time, event, IPTW=1) {
	data.frame(est_effect=est_effect, treatment=treatment, time=time, event=event, IPTW=IPTW) %>%
		mutate(do_treat = est_effect < 0) %>% # treat if the log-relative risk is negative
		filter(do_treat == treatment) %>% # the "lucky" individuals
		arrange(-time) %>%
		mutate(weighted_n_alive = cumsum(IPTW)) %>%
		arrange(time) %>%
		mutate(km_surv = cumprod(1-(event*IPTW)/weighted_n_alive)) %>%
		mutate(dt = lead(time) - time) %>%
		mutate(rest_mean_survival = km_surv*dt) %>%
		pull(rest_mean_survival) %>% 
		sum(na.rm=T) # an estimate of the restricted mean survival time of the lucky patients
}

true_value = function(est_effect, treated_mean, control_mean) {
	do_treat = est_effect < 0 # treat if the log-relative risk is negative
	obs_mean = do_treat*treated_mean + (1-do_treat)*control_mean
	mean(obs_mean)
}

random_metric = function(){ return(0)}