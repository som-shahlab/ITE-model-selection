library(Matching)
library(magrittr)
library(tidyverse)

# Refactor this into methods that estimate \check\tau and then losses? (except for value)

### MSE Methods ###

mse = function(estimates, obs_outcome) {
	mean((estimates-obs_outcome)^2) # mean squared error
}

transformed_mse = function(est_effect, obs_treatment, obs_outcome, weights=1) {
	# transformed_outcome = mean((weights*obs_outcome*(2*obs_treatment - 1) - est_effect)^2)
	transformed_outcome = mse(est_effect, weights*obs_outcome*(2*obs_treatment - 1))
}

matching_mse = function(est_effect, obs_treatment, obs_outcome, subject, match) { # need to do for all test sets
	# for each pair, calculate the difference in outcomes (in the right direction)
	data.frame(subject=subject, match=match, treatment=obs_treatment, 
			   subject_outcome=obs_outcome, effect=est_effect) %>%
		inner_join(data.frame(match=subject, match_outcome=obs_outcome), by="match") %>%
		mutate(delta = (2*treatment-1)*(subject_outcome - match_outcome)) %>%
		mutate(sq_error = (delta-est_effect)^2) %>%
		pull(sq_error) %>%
		mean()
}

# This one should not work
est_te_strata_mse = function(est_effect, obs_treatment, obs_outcome, n_strata=10) {
	data.frame(est_effect, obs_outcome, obs_treatment) %>%
		mutate(strata = ntile(est_effect, n_strata)) %>%
		dplyr::group_by(strata) %>%
		dplyr::summarize(mean_hte_estimate = mean(est_effect),
				  mean_treated_outcome = sum(obs_outcome*obs_treatment)/sum(obs_treatment),
				  mean_control_outcome = sum(obs_outcome*!obs_treatment)/sum(!obs_treatment),
				  n_in_strata = n()) %>%
		dplyr::mutate(te_estimate = mean_treated_outcome - mean_control_outcome,
			   	  	  error = n_in_strata*(mean_hte_estimate - te_estimate)^2) %>%
		filter(!is.na(error)) %>%
		pull(error) %>%
		mean()
}

### Value Methods ###

true_hte_value = function(est_effect, true_effect, true_mean, cutoff=0) {
	do_treat = est_effect >= cutoff
	mean(true_mean + true_effect*(2*do_treat - 1) / 2)
}

# c_benefit_k = function(est_effect, obs_treatment, obs_outcome, cutoff=0, weights=1) {
# 	do_treat = est_effect >= cutoff
# 	weighted_outcome = weights*obs_outcome
# 	(mean(weighted_outcome[do_treat & obs_treatment]) - mean(weighted_outcome[do_treat & !obs_treatment]))*sum(do_treat)/length(est_effect) -
# 	(mean(weighted_outcome[!do_treat & obs_treatment]) - mean(weighted_outcome[!do_treat & !obs_treatment]))*sum(!do_treat)/length(est_effect)
# }

gain = function(est_effect, obs_treatment, obs_outcome, cutoff=0, weights=1) {
	do_treat = est_effect >= cutoff
	weighted_outcome = weights*obs_outcome
	(mean(weighted_outcome[do_treat & obs_treatment]) - mean(weighted_outcome[do_treat & !obs_treatment]))*sum(do_treat)/length(est_effect)
}

value = function(est_effect, obs_treatment, obs_outcome, cutoff=0, weights=1) {
	do_treat = est_effect >= cutoff
	weighted_outcome = weights*obs_outcome
	sum(weighted_outcome[do_treat==obs_treatment])/(length(est_effect))
}

# value_max = function(est_effect, obs_treatment, obs_outcome, weights=1) {
# 	weighted_outcome = weights*obs_outcome
# 	data.frame(est_effect, obs_treatment, weighted_outcome) %>%
# 		arrange(-est_effect) %>%
# 		mutate(Yt_lucky = cumsum(weighted_outcome*obs_treatment)) %>%
# 		arrange(est_effect) %>%
# 		mutate(Yc_lucky = cumsum(weighted_outcome*!obs_treatment)) %>%
# 		mutate(value = (Yt_lucky + Yc_lucky)/n()) %>%
# 		filter(!is.nan(value)) %>%
# 		pull(value) %>% max
# }

### AUC-type Methods ###

c_benefit = function(est_effect, obs_treatment, obs_outcome) {
	match = Match(Tr=obs_treatment, 
				  X=est_effect,
				  replace=T, estimand="ATT")
	delta = obs_outcome[match$index.treated] - obs_outcome[match$index.control]
	ranked_pairs = data.frame(delta, est_effect=est_effect[match$index.treated]) %>%
		arrange(est_effect) %>% # from smallest to biggest 
		mutate(effect_rank=row_number(), trash_join_var=1)
	inner_join(ranked_pairs, ranked_pairs, by="trash_join_var") %>%
		mutate(concordant=((effect_rank.x<effect_rank.y) & (delta.x<delta.y))) %>%
		pull(concordant) %>% 
		mean()
}

uplift = function(est_effect, obs_treatment, obs_outcome, weights=1) {
	weighted_outcome = weights*obs_outcome
	data.frame(est_effect, obs_treatment, weighted_outcome) %>%
		arrange(-est_effect) %>%
		mutate(Yt = cumsum(weighted_outcome*obs_treatment),
			   Yc = cumsum(weighted_outcome*(!obs_treatment)),
			   Nt = cumsum(obs_treatment),
			   Nc = cumsum(!obs_treatment)) %>%
		mutate(uplift = (Yt/Nt - Yc/Nc)*(Nt+Nc)/length(est_effect)) %>%
		filter(!is.nan(uplift)) %>%
		sum()  
}