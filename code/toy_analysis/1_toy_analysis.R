

#=============================================================================#
# setting up workspace
#=============================================================================#
rm(list = ls())
setwd('.')

#=============================================================================#
# creating toy data sets
#=============================================================================#

pathogen_a = seq(0, 1000, by = 1)
pathogen_b1 = 100
pathogen_b2 = 1000
pathogen_b3 = 10000

p_prime_a1 = pathogen_a / (pathogen_a+pathogen_b1)
p_prime_a2 = pathogen_a / (pathogen_a+pathogen_b2)
p_prime_a3 = pathogen_a / (pathogen_a+pathogen_b3)
plot(p_prime_a1, ylim = c(0,1), col = 'deeppink3', type = 'l', lwd = 3)
lines(p_prime_a2, col = 'mediumseagreen', type = 'l', lwd = 3)
lines(p_prime_a3, col = 'slateblue', type = 'l', lwd = 3)


## braga values
# se = c(0.813, 1, 0.583, 0.809, 0.756, 0.286)
# se[which(se == 1)] = 0.999
# sp = c(0.109, 0.014, 0.519, 0.580, 0.635, 0.973)
# se_sp_mat = matrix(NA, ncol = 2, nrow = length(se))
# se_sp_mat[,1] = se
# se_sp_mat[,2] = sp
# 
# se_sp_mat = se_sp_mat[order(se),]

## random values
# 6x6
# se = c(0.2, 0.4, 0.6, 0.8)
# sp = c(0.2, 0.4, 0.6, 0.8)

# 5x5
se = c(0.8, 0.5, 0.2)
sp = c(0.2, 0.5, 0.8)


# se = c(0.1, 0.5, 0.9)
# sp = c(0.1, 0.5, 0.9)

# se = seq(0.1, 0.9, by = 0.3)
se = c(0.99, se)
# sp = seq(0.1, 0.9, by = 0.3)
sp = c(sp, 0.99)
se_sp_mat = expand.grid(se, sp)
# se_sp_mat[100,] = c(0.286, 0.973)
names(se_sp_mat) = c('se', 'sp')
se_sp_mat = as.matrix(se_sp_mat, ncol = 2)


# se_sp_mat = matrix(c(0.286, 0.973), nrow = 1)

#=============================================================================#
# functions
#=============================================================================#

# se and sp constraints 1 & 2 --> 1 = se > p'; 2 = se + sp cannot equal 1
constrain_diag_1_2 = function(p_prime_in, diag_mat_in, browse = F){
  if(browse) browser()
  return(which((diag_mat_in[,1] > p_prime_in) & diag_mat_in[,1] + diag_mat_in[,2] != 1))
}


# se and sp constraint 3 --> numerator and denomiator of p must either both be 
# positive or negative to result in a postive p value
constrain_diag_3 = function(p_prime_in, diag_mat_in, browse = F){
  if(browse) browser()
  
  if(length(diag_mat_in) == 0){
    diag_mat_in = matrix(NA, ncol= 2)
  }
  
  if(length(nrow(diag_mat_in)) == 0){
    diag_mat_in = matrix(diag_mat_in, ncol = 2)
  }
  
  a = sapply(1:nrow(diag_mat_in), function(ff) p_prime_in - 1 + diag_mat_in[ff,2])
  b = sapply(1:nrow(diag_mat_in), function(ff) diag_mat_in[ff,1] - 1 + diag_mat_in[ff,2])
  
  # out_return = which(a/b <= 1 & a/b >= 0)
  # if(sum(out_return) == 0){
  #   out_return = NA
  # }

  # return(out_return)
  return(which(a/b <= 1 & a/b >= 0))
}


# calculate pZ function 
calc_pZ = function(se_in, sp_in, p_prime_in, browse = F){
  if(browse) browser()
  return(
    sapply(1:length(se_in), function(ff) 
      ((p_prime_in - 1 + sp_in[ff]) / (sp_in[ff] - 1 + se_in[ff])))
  )  
}


# function to calculate the number of missed Zika infections
# if output is positive --> there were Zika cases missed
# if output is negative --> there were dengue/chik cases missed
calc_missed_zika = function(c_Z_in, c_O_in, p_z_in){
  cases_missed_out = ((c_Z_in + c_O_in) * p_z_in) - c_Z_in
  return(cases_missed_out)
}


#=============================================================================#
# calculating allowable se and sp and estimating p_Z
#=============================================================================#

p_a1 = list()
length(p_a1) = length(p_prime_a1)
p_a2 = list()
length(p_a2) = length(p_prime_a2)
p_a3 = list()
length(p_a3) = length(p_prime_a3)

  
for(ii in 1:length(p_prime_a2)){
  diag_1_2 = constrain_diag_1_2(p_prime_a1[ii], se_sp_mat)
  diag_3 = matrix(se_sp_mat[diag_1_2,], ncol = 2)[constrain_diag_3(p_prime_a1[ii], se_sp_mat[diag_1_2,]),]

  p_a1[[ii]] = list(calc_pZ(se_in = matrix(diag_3, ncol = 2)[,1],
                            sp_in = matrix(diag_3, ncol = 2)[,2],
                            p_prime_in = p_prime_a1[ii], browse = F), diag_3)
  
  diag_1_2 = constrain_diag_1_2(p_prime_a2[ii], se_sp_mat)
  diag_3 = matrix(se_sp_mat[diag_1_2,], ncol = 2)[constrain_diag_3(p_prime_a2[ii], se_sp_mat[diag_1_2,]),]
  
  p_a2[[ii]] = list(calc_pZ(se_in = matrix(diag_3, ncol = 2)[,1], 
                            sp_in = matrix(diag_3, ncol = 2)[,2], 
                            p_prime_in = p_prime_a2[ii]), diag_3)
  
  diag_1_2 = constrain_diag_1_2(p_prime_a3[ii], se_sp_mat)
  diag_3 = matrix(se_sp_mat[diag_1_2,], ncol = 2)[constrain_diag_3(p_prime_a3[ii], se_sp_mat[diag_1_2,]),]

  p_a3[[ii]] = list(calc_pZ(se_in = matrix(diag_3, ncol = 2)[,1],
                            sp_in = matrix(diag_3, ncol = 2)[,2],
                            p_prime_in = p_prime_a3[ii]), diag_3)
  
}


#=============================================================================#
# plotting observed versus re-estimated
#=============================================================================#

se_sp_mat

df_list1 = list()
length(df_list1)= nrow(se_sp_mat)

df_list2 = list()
length(df_list2)= nrow(se_sp_mat)

df_list3 = list()
length(df_list3)= nrow(se_sp_mat)


for(ii in 1:length(df_list1)){
  df_list1[[ii]] = rep(NA, length(p_a2))
  df_list2[[ii]] = rep(NA, length(p_a2))
  df_list3[[ii]] = rep(NA, length(p_a2))
}



for(ii in 1:length(p_a2)){
  
  #b1
  ind_tmp = which(se_sp_mat[,1] %in% matrix(p_a1[[ii]][[2]], ncol = 2)[,1] & 
                    se_sp_mat[,2] %in% matrix(p_a1[[ii]][[2]], ncol = 2)[,2])
  if(sum(ind_tmp != 0)){
    for(tt in 1:length(ind_tmp)){
      df_list1[[ind_tmp[tt]]][ii] = calc_missed_zika(c_Z_in = pathogen_a[ii], 
                                                    c_O_in = pathogen_b1, 
                                                    p_z_in = p_a1[[ii]][[1]][tt])
    }
  }
  
  # b2
  ind_tmp = which(se_sp_mat[,1] %in% matrix(p_a2[[ii]][[2]], ncol = 2)[,1] & 
                    se_sp_mat[,2] %in% matrix(p_a2[[ii]][[2]], ncol = 2)[,2])
  if(sum(ind_tmp != 0)){
    for(tt in 1:length(ind_tmp)){
      df_list2[[ind_tmp[tt]]][ii] = calc_missed_zika(c_Z_in = pathogen_a[ii], 
                                                    c_O_in = pathogen_b2, 
                                                    p_z_in = p_a2[[ii]][[1]][tt])
    }
  }
  
  #b3
  ind_tmp = which(se_sp_mat[,1] %in% matrix(p_a3[[ii]][[2]], ncol = 2)[,1] & 
                    se_sp_mat[,2] %in% matrix(p_a3[[ii]][[2]], ncol = 2)[,2])
  if(sum(ind_tmp != 0)){
    for(tt in 1:length(ind_tmp)){
      df_list3[[ind_tmp[tt]]][ii] = calc_missed_zika(c_Z_in = pathogen_a[ii], 
                                                    c_O_in = pathogen_b3, 
                                                    p_z_in = p_a3[[ii]][[1]][tt])
    }
  }
}


revised_a1 = df_list1
revised_a2 = df_list2
revised_a3 = df_list3
for(ii in 1:length(df_list1)){
  revised_a1[[ii]] = df_list1[[ii]] + pathogen_a
  revised_a2[[ii]] = df_list2[[ii]] + pathogen_a
  revised_a3[[ii]] = df_list3[[ii]] + pathogen_a
}


pdf('../output/results_final/toy_analysis_revised_estimates_se_sp_set.pdf',
    height = 5.5, width = 6.5)
# layout(matrix(1:36, 6, 6))
# layout(matrix(1:25, 5, 5))
layout(matrix(1:16, 4, 4))
# layout(matrix(1:9, 3, 3))
# layout(matrix(1:6, 2, 3))
par(mar = c(0.25, 0.25, 0.25, 0.25), oma = c(5,5,2,6))

for(ii in 1:length(revised_a1)){
  
  mat_tmp = matrix(NA, nrow = length(pathogen_a), ncol = 4)
  mat_tmp[,1] = pathogen_a
  mat_tmp[,2] = revised_a1[[ii]]
  mat_tmp[,3] = revised_a2[[ii]]
  mat_tmp[,4] = revised_a3[[ii]]
  
  # if(ii %in% c(6, 12, 18, 24, 30, 36)){
  if(ii %in% c(4, 8, 12, 16)){  
    if(ii == 4){
      plot(-100, -100, xlim = c(0, 1000), ylim = c(0,5000), las = 1, xaxt = 'n')
      axis(side = 1, at = c(0, 500, 1000), labels = c(0, 500, 1000), las = 2)
    }else{
      plot(-100, -100, xlim = c(0, 1000), ylim = c(0,5000), yaxt = 'n', xaxt = 'n')
      axis(side = 1, at = c(0, 500, 1000), labels = c(0, 500, 1000), las = 2)
    }
  }else if(ii %in% c(1:4)){
    plot(-100, -100, xlim = c(0, 1000), ylim = c(0,5000), xaxt = 'n', las = 1)
  }else{
    plot(-100, -100, xlim = c(0, 1000), ylim = c(0,5000), xaxt = 'n', yaxt = 'n')
  }
  
  if(ii == 1){
    legend('topleft', legend = c('Disease B = 100', 
                                          'Disease B = 1000',
                                          'Disease B = 10000',
                                 '1:1 line'),
           lty = c(1,1,1,3), lwd = 2, col = c('deeppink3', 'mediumseagreen', 'slateblue', rgb(0,0,0,0.5)), bty = 'n', cex = 0.85)
    
    mtext(side = 3, text = expression(italic('sp = ') * 0.2))
  }
  
  if(ii == 5){
    mtext(side = 3, text = expression(italic('sp = ') * 0.5))
  }
  
  if(ii == 9){
    mtext(side = 3, text = expression(italic('sp = ') * 0.8))
  }
  
  if(ii == 13){
    mtext(side = 3, text = expression(italic('sp = ') * 0.99))
    mtext(side = 4, text = expression(italic('se = ') * 0.99), las = 1, line = 0.25)
  }
  
  if(ii == 14){
    mtext(side = 4, text = expression(italic('se = ') * 0.8), las = 1, line = 0.25)
  }
  
  if(ii == 15){
    mtext(side = 4, text = expression(italic('se = ') * 0.5), las = 1, line = 0.25)
  }
  
  if(ii == 16){
    mtext(side = 4, text = expression(italic('se = ') * 0.2), las = 1, line = 0.25)
  }
  
  abline(1,1, col = rgb(0,0,0,0.5), lwd = 2, lty = 3)
  
  lines(mat_tmp[,1], mat_tmp[,2], col = 'deeppink3', lwd = 2)
  # points(mat_tmp[,1], mat_tmp[,2], pch = 16, col = 'deeppink3', cex = 1.2)
  
  lines(mat_tmp[,1], mat_tmp[,3], col = 'mediumseagreen', lwd = 2)
  # points(mat_tmp[,1], mat_tmp[,3], pch = 16, col = 'mediumseagreen', cex = 1.2)
  
  lines(mat_tmp[,1], mat_tmp[,4], col = 'slateblue', lwd = 2)
  # points(mat_tmp[,1], mat_tmp[,4], pch = 16, col = 'slateblue', cex = 1.2)
  
  # mtext(side = 3, text = paste0('se=', se_sp_mat[ii,1], ', sp=', se_sp_mat[ii,2]), cex = 0.9)
}
mtext(side =1, outer = T, 'Reported incidence of disease A', line = 3.5)
mtext(side =2, outer = T, 'Revised incidence of disease A', line = 3.5)

dev.off()

