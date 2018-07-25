
# reverse the order of a df (sort without reordering)
reverse_rows = function(df) df[dim(df)[1]:1,]

match_tau = function(w, y, matches) (2*w-1)*(y-y[matches])
trans_tau = function(w, y, iptw) iptw*y*(2*w-1)

# "mean decision cost"
mdc = function(tau, tau_hat) -mean(tau*(tau_hat > 0))

####################################################
######################## μ-risk ####################
####################################################

#' @export
mse = function(y, y_hat) mean((y_hat-y)^2)
#' @export
bundle_mse = function(bundle, y_hat, tau_hat) bundle %$% mse(y, y_hat)

#' @export
wmse = function(y, iptw, y_hat) sum(iptw*(y_hat-y)^2)/sum(iptw)
#' @export
bundle_wmse = function(bundle, y_hat, tau_hat) bundle %$% wmse(y, iptw, y_hat)



####################################################
######################## τ-risk ####################
####################################################

#' @export
dr_mse = function(w, y, p,  mu0, mu1, tau_hat) sum(((y - (w-p)*tau_hat)/(p*(1-p)))^2)
#' @export
bundle_dr_mse = function(bundle, y_hat, tau_hat) bundle %$% dr_mse(w, y, p,  mu0, mu1, tau_hat)

#' @export
r_objective = function(weight, pseudo_outcome, tau_hat) wmse(pseudo_outcome, weight, tau_hat)
#' @export
bundle_r_objective = function(bundle, y_hat, tau_hat) bundle %$% r_objective(weight, pseudo_outcome, tau_hat)

#' @export
match_mse = function(matches, w, y, tau_hat) match_tau(w,y,matches) %>% mse(tau_hat)
#' @export
bundle_match_mse = function(bundle, y_hat, tau_hat) bundle %$% match_mse(matches, w, y, tau_hat)

#' @export
trans_mse = function(w, y, iptw, tau_hat) trans_tau(w,y,iptw) %>% mse(tau_hat)
#' @export
bundle_trans_mse = function(bundle, y_hat, tau_hat) bundle %$% trans_mse(w, y, iptw, tau_hat)

#' @export
match_mdc = function(matches, w, y, tau_hat) match_tau(w,y,matches) %>% mdc(tau_hat)
#' @export
bundle_match_mdc = function(bundle, y_hat, tau_hat) bundle %$% match_mdc(matches, w, y, tau_hat)

# generalized gain
#' @export
trans_mdc =  function(w, y, iptw, tau_hat) trans_tau(w,y,iptw) %>% mdc(tau_hat)
#' @export
bundle_trans_mdc = function(bundle, y_hat, tau_hat) bundle %$% trans_mdc(w, y, iptw, tau_hat)

#####################################################
######################### value #####################
#####################################################

#' @export
ip_value = function(w, y, iptw, tau_hat) sum((iptw*y)[(tau_hat > 0)==w])/length(w)
#' @export
bundle_ip_value = function(bundle, y_hat, tau_hat) bundle %$% ip_value(w, y, iptw, tau_hat)

#' @export
dml_value = function(w, y, iptw, mu0, mu1, tau_hat) {
	gamma = mu1 - mu0 + trans_tau(w, y, iptw) - w*trans_tau(w, mu1, iptw) + (1-w)*trans_tau(w, mu0, iptw)
	mean((2*(tau_hat>0)-1)*gamma)
}
#' @export
bundle_dml_value = function(bundle, y_hat, tau_hat) bundle %$% dml_value(w, y, iptw, mu0, mu1, tau_hat)

# old gain
#' @export
gain = function(w, y, tau_hat) {
	d = tau_hat>0 
	(mean(y[d & w]) - mean(y[d & !w]))*sum(d)/length(w)}
#' @export
bundle_gain = function(bundle, y_hat, tau_hat) bundle %$% gain(w, y, tau_hat)

#####################################################
######################### AUC #######################
#####################################################

#' @export
c_benefit = function(w, y, tau_hat) {
	match = Match(Tr=w, 
				  X=tau_hat,
				  replace=T, estimand="ATT")
	delta = y[match$index.treated] - y[match$index.control]
	ranked_pairs = data.frame(delta, tau_hat=tau_hat[match$index.treated]) %>%
		arrange(tau_hat) %>% # from smallest to biggest 
		mutate(effect_rank=row_number(), trash_join_var=1)
	inner_join(ranked_pairs, ranked_pairs, by="trash_join_var") %>%
		mutate(concordant=((effect_rank.x<effect_rank.y) & (delta.x<delta.y))) %>%
		pull(concordant) %>% 
		mean()
}
#' @export
bundle_c_benefit = function(bundle, y_hat, tau_hat) bundle %$% c_benefit(w, y, tau_hat)

#' @export
qini = function(w, y, iptw, tau_hat) {
	weighted_y = iptw*y
	data.frame(tau_hat, w, y) %>%
		arrange(-tau_hat) %>%
		mutate(Yt = cumsum(y*w),
			   Yc = cumsum(y*(!w)),
			   Nt = cumsum(w),
			   Nc = cumsum(!w)) %>%
		mutate(uplift = (Yt/Nt - Yc/Nc)*(Nt+Nc)/n()) %>%
		filter(!is.nan(uplift)) %>%
		sum()  
}
#' @export
bundle_qini = function(bundle, y_hat, tau_hat) bundle %$% qini(w, y, iptw, tau_hat)

# ip_value_auc = function(w, y, iptw, tau_hat) {
# 	data.frame(tau_hat, w, y, iptw) %>%
# 		arrange(-tau_hat) %>%
# 		mutate(Yt_lucky = cumsum(iptw*y*w)) %>%
# 		reverse_rows() %>% # arrange(tau_hat), but faster
# 		mutate(Yc_lucky = cumsum(iptw*y*!w)) %>%
# 		mutate(value = (Yt_lucky + Yc_lucky)/n()) %>%
# 		filter(!is.nan(value)) %>%
# 		pull(value) %>% 
# 		sum()
# }
# ip_bundle_value_auc = function(bundle, y_hat, tau_hat) bundle %$% ip_value_auc(w, y, iptw, tau_hat)

#####################################################
######################### Random #######################
#####################################################

#' @export
random_metric = function() return(0)
#' @export
bundle_random = function(bundle, y_hat, tau_hat) random_metric()