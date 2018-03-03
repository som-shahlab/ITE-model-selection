library(Matching)
library(magrittr)
library(tidyverse)

####################################################
############# FRAMEWORK METHODS ####################
####################################################

### Estimators for \check\tau
est_effect_covariate_matching = function(treatment, outcome, subject, match) {
	data.frame(subject=subject, match=match, treatment=treatment, 
			   subject_outcome=outcome) %>%
		inner_join(data.frame(match=subject, match_outcome=outcome), by="match") %>%
		mutate(est_effect_test = (2*treatment-1)*(subject_outcome - match_outcome)) %>%
		pull(est_effect_test)
}

est_effect_transformed_outcome = function(treatment, outcome, weights) {
	weights*outcome*(2*treatment - 1)
}

### Loss functions
loss_squared_error = function(truth, estimate) {
	mean((estimate-truth)^2)
}

loss_decision = function(truth, estimate, cutoff=0) {
	mean(-truth*(estimate>=cutoff))
}

#####################################################
############# NON-FRAMEWORK METHODS #################
#####################################################

# Rob's broken method
est_te_strata = function(est_effect, treatment, outcome, n_strata=10) {
	data.frame(est_effect, outcome, treatment) %>%
		mutate(strata = ntile(est_effect, n_strata)) %>%
		dplyr::group_by(strata) %>%
		dplyr::summarize(mean_hte_estimate = mean(est_effect),
				  mean_treated_outcome = sum(outcome*treatment)/sum(treatment),
				  mean_control_outcome = sum(outcome*!treatment)/sum(!treatment),
				  n_in_strata = n()) %>%
		dplyr::mutate(te_estimate = mean_treated_outcome - mean_control_outcome,
			   	  	  error = n_in_strata*(mean_hte_estimate - te_estimate)^2) %>%
		filter(!is.na(error)) %>%
		pull(error) %>%
		mean()
}

### Value Methods ###
true_value = function(est_effect, true_effect, true_mean, cutoff=0) {
	do_treat = est_effect >= cutoff
	mean(true_mean + true_effect*(2*do_treat - 1) / 2)
}

value = function(est_effect, treatment, outcome, cutoff=0, weights=1) {
	do_treat = est_effect >= cutoff
	weighted_outcome = weights*outcome
	sum(weighted_outcome[do_treat==treatment])/(length(est_effect))
}

gain = function(est_effect, treatment, outcome, cutoff=0, weights=1) {
	do_treat = est_effect >= cutoff
	weighted_outcome = weights*outcome
	(mean(weighted_outcome[do_treat & treatment]) - mean(weighted_outcome[do_treat & !treatment]))*sum(do_treat)/length(est_effect)
}

# c_benefit_k = function(est_effect, treatment, outcome, cutoff=0, weights=1) {
# 	do_treat = est_effect >= cutoff
# 	weighted_outcome = weights*outcome
# 	(mean(weighted_outcome[do_treat & treatment]) - mean(weighted_outcome[do_treat & !treatment]))*sum(do_treat)/length(est_effect) -
# 	(mean(weighted_outcome[!do_treat & treatment]) - mean(weighted_outcome[!do_treat & !treatment]))*sum(!do_treat)/length(est_effect)
# }

# value_max = function(est_effect, treatment, outcome, weights=1) {
# 	weighted_outcome = weights*outcome
# 	data.frame(est_effect, treatment, weighted_outcome) %>%
# 		arrange(-est_effect) %>%
# 		mutate(Yt_lucky = cumsum(weighted_outcome*treatment)) %>%
# 		arrange(est_effect) %>%
# 		mutate(Yc_lucky = cumsum(weighted_outcome*!treatment)) %>%
# 		mutate(value = (Yt_lucky + Yc_lucky)/n()) %>%
# 		filter(!is.nan(value)) %>%
# 		pull(value) %>% max
# }

### AUC-type Methods ###

c_benefit = function(est_effect, treatment, outcome) {
	match = Match(Tr=treatment, 
				  X=est_effect,
				  replace=T, estimand="ATT")
	delta = outcome[match$index.treated] - outcome[match$index.control]
	ranked_pairs = data.frame(delta, est_effect=est_effect[match$index.treated]) %>%
		arrange(est_effect) %>% # from smallest to biggest 
		mutate(effect_rank=row_number(), trash_join_var=1)
	inner_join(ranked_pairs, ranked_pairs, by="trash_join_var") %>%
		mutate(concordant=((effect_rank.x<effect_rank.y) & (delta.x<delta.y))) %>%
		pull(concordant) %>% 
		mean()
}

qini = function(est_effect, treatment, outcome, weights=1) {
	weighted_outcome = weights*outcome
	data.frame(est_effect, treatment, weighted_outcome) %>%
		arrange(-est_effect) %>%
		mutate(Yt = cumsum(weighted_outcome*treatment),
			   Yc = cumsum(weighted_outcome*(!treatment)),
			   Nt = cumsum(treatment),
			   Nc = cumsum(!treatment)) %>%
		mutate(uplift = (Yt/Nt - Yc/Nc)*(Nt+Nc)/n()) %>%
		filter(!is.nan(uplift)) %>%
		sum()  
}

value_auc = function(est_effect, treatment, outcome, weights=1) {
	weighted_outcome = weights*outcome
	data.frame(est_effect, treatment, weighted_outcome) %>%
		arrange(-est_effect) %>%
		mutate(Yt_lucky = cumsum(weighted_outcome*treatment)) %>%
		arrange(est_effect) %>%
		mutate(Yc_lucky = cumsum(weighted_outcome*!treatment)) %>%
		mutate(value = (Yt_lucky + Yc_lucky)/n()) %>%
		filter(!is.nan(value)) %>%
		pull(value) %>% 
		sum()
}