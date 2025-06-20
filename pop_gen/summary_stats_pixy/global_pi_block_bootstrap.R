#library(dplyr)

boot_pi = function(diffs_vec = NULL, comps_vec = NULL, n_boot = 1000){

    if( length(diffs_vec) != length(comps_vec) ){
        stop("Error: vectors are of different lengths")
    }
    
    if( sum(is.na(diffs_vec) != is.na(comps_vec)) > 0){
        stop("Error: vectors have unmatched NA values")
    }
    
    
    last_ind=length(diffs_vec)
    pi_vec = vector(mode = "numeric", length = n_boot)
    
    for(i in 1:n_boot){
        boot_inds = sample(1:last_ind, replace = T)
        
        pi_vec[i] = sum(diffs_vec[boot_inds], na.rm = T)/sum(comps_vec[boot_inds], na.rm = T)
    }
    return(pi_vec)
}

pi_dat = read.table("data/Nf/pixy/window_10kb_single_pop/pixy_pi.txt", header = T)
head(pi_dat)

pi_boots = boot_pi(diffs_vec = pi_dat$count_diffs, comps_vec = pi_dat$count_comparisons, n_boot = 10000)
#pi_boots
hist(pi_boots)

pi_quants = quantile(pi_boots, probs = c(0.025, 0.5, 0.975))
pi_quants
#        2.5%         50%       97.5% 
# 0.001636438 0.001784834 0.001938701 
mean(pi_boots)-pi_quants[1]
#0.0001497085
pi_quants[3]-mean(pi_boots)
#0.0001525545


# Nd
pi_dat = read.table("data/Nd/pixy/window_10kb_single_pop/pixy_pi.txt", header = T)

pi_boots = boot_pi(diffs_vec = pi_dat$count_diffs, comps_vec = pi_dat$count_comparisons, n_boot = 10000)
#pi_boots
hist(pi_boots)

pi_quants = quantile(pi_boots, probs = c(0.025, 0.5, 0.975))
pi_quants
#        2.5%         50%       97.5% 
# 0.003769640 0.003866358 0.003967494 
mean(pi_boots)-pi_quants[1]
#0.00009719266
pi_quants[3]-mean(pi_boots)
#0.0001006623


# Nc
pi_dat = read.table("data/Nc/pixy/window_10kb_single_pop/pixy_pi.txt", header = T)

pi_boots = boot_pi(diffs_vec = pi_dat$count_diffs, comps_vec = pi_dat$count_comparisons, n_boot = 10000)
#pi_boots
hist(pi_boots)

pi_quants = quantile(pi_boots, probs = c(0.025, 0.5, 0.975))
pi_quants
#       2.5%         50%       97.5% 
# 0.003361187 0.003473057 0.003592182 
mean(pi_boots)-pi_quants[1]
#0.0001124276
pi_quants[3]-mean(pi_boots)
#0.0001185673

pi_boots[2]-pi_quants[1]
#0.0001179435
pi_quants[3]-pi_boots[2]
#0.0001130513

