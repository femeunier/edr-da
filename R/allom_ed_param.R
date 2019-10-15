#' Get All ED2 allometry parameters
#'
#' @export
#' 

allom_ed_param <- function(allom=3){
  
  npft = 17
  
  odead_small = c(-1.1138270, 2.4404830, 2.1806320 )
  odead_large = c( 0.1362546, 2.4217390, 6.9483532 )
  ndead_small = c(-1.2639530, 2.4323610, 1.8018010 )
  ndead_large = c(-0.8346805, 2.4255736, 2.6822805 )
  nleaf       = c( 0.0192512, 0.9749494, 2.5858509 )
  ncrown_area = c( 0.1184295, 1.0521197            )
  
  pft = list()
  
  pft$is_tropical = logical(npft)
  pft$is_tropical[1:4]   = TRUE
  pft$is_tropical[5:13]  = FALSE
  pft$is_tropical[14:17] = TRUE
  
  pft$is_liana = integer(npft)
  pft$is_liana[1:16]  = FALSE
  pft$is_liana[17]    = TRUE
  
  pft$is_grass = integer(npft)
  pft$is_grass[1]     = TRUE
  pft$is_grass[2:4]   = FALSE
  pft$is_grass[5]     = TRUE
  pft$is_grass[6:11]  = FALSE
  pft$is_grass[12:16] = TRUE
  pft$is_grass[17]    = FALSE
  
  pft$rho = numeric(npft)
  pft$rho[1]     = 0.20   
  pft$rho[2]     = 0.53   
  pft$rho[3]     = 0.71  
  pft$rho[4]     = 0.90   
  pft$rho[5]     = 0.20   
  pft$rho[6]     = 0.45  
  pft$rho[7]     = 0.57                
  pft$rho[8]     = 0.45   
  pft$rho[9]     = 0.42  
  pft$rho[10]    = 0.74   
  pft$rho[11]    = 0.70   
  pft$rho[12:16] = 0.20
  pft$rho[17]    = 0.46   
  
  pft$b1Ht = numeric(npft)
  pft$b1Ht[1:4]   = 0.0
  pft$b1Ht[5]     = 0.4778
  pft$b1Ht[6]     = 27.14
  pft$b1Ht[7]     = 27.14
  pft$b1Ht[8]     = 22.79
  pft$b1Ht[9]     = 22.6799
  pft$b1Ht[10]    = 25.18
  pft$b1Ht[11]    = 23.3874
  pft$b1Ht[12:13] = 0.4778
  pft$b1Ht[14:15] = 0.0
  pft$b1Ht[16]    = 0.0
  pft$b1Ht[17]    = 0.0
  
  pft$b2Ht = numeric(npft)
  pft$b2Ht[1:4]   =  0.00
  pft$b2Ht[5]     = -0.75
  pft$b2Ht[6]     = -0.03884
  pft$b2Ht[7]     = -0.03884
  pft$b2Ht[8]     = -0.04445
  pft$b2Ht[9]     = -0.06534
  pft$b2Ht[10]    = -0.04964
  pft$b2Ht[11]    = -0.05404
  pft$b2Ht[12:13] = -0.75
  pft$b2Ht[14:15] =  0.00
  pft$b2Ht[16]    =  0.00
  pft$b2Ht[17]    =  0.00
  
  pft$hgt_ref = numeric(npft)
  pft$hgt_ref[1:5]   = 0.0
  pft$hgt_ref[6:11]  = 1.3
  pft$hgt_ref[12:15] = 0.0
  pft$hgt_ref[16]    = 0.0
  pft$hgt_ref[17]    = 0.0
  
  if (allom<2) {
    for(ipft in seq(1,npft)){
      if (pft$is_tropical[ipft]){
        pft$b1Ht   [ipft] = 0.37 * log(10.0)
        pft$b2Ht   [ipft] = 0.64
      }
    }
    
  } else {
    for(ipft in seq(1,npft)){
      if (pft$is_tropical[ipft]){
        pft$b1Ht   [ipft] = 0.0352
        pft$b2Ht   [ipft] = 0.694
        pft$hgt_ref[ipft] = 61.7
      }
    }
  }
  
  pft$b1Ht[17] = 0.1136442
  pft$b2Ht[17] = 0.8675

  pft$b1Ca = numeric(npft)
  pft$b1Ca[1:4]   = exp(-1.853) * exp(pft$b1Ht[1:4]) ** 1.888
  pft$b1Ca[5:13]= 2.490154
  pft$b1Ca[14:17] = exp(-1.853) * exp(pft$b1Ht[14:17]) ** 1.888
  
  pft$b2Ca = numeric(npft)
  pft$b2Ca[1:4]   = pft$b2Ht[1:4] * 1.888
  pft$b2Ca[5:13]  = 0.8068806
  pft$b2Ca[14:17] = pft$b2Ht[14:17] * 1.888
  
  if (allom < 2) {
    for(ipft in seq(1,npft)){
      if (pft$is_tropical[ipft]){
        pft$b1Ca[ipft] = exp(-1.853) * exp(pft$b1Ht[ipft]) ** 1.888
        pft$b2Ca[ipft] = pft$b2Ht[ipft] * 1.888
      }
    }
  } else {
    for(ipft in seq(1,npft)){
      if (pft$is_tropical[ipft]){
        pft$b1Ca[ipft] = exp(ncrown_area[1])
        pft$b2Ca[ipft] = ncrown_area[2]
      }
    }
  }
  
  pft$b2Ca[17] = 1.26254364 
  
  clump_grass = 1
  clump_tree = 0.8
  
  pft$clumping_factor[1]     = clump_grass
  pft$clumping_factor[2:4]   = clump_tree 
  pft$clumping_factor[5]     = 8.400e-1
  pft$clumping_factor[6:8]   = 7.350e-1
  pft$clumping_factor[9:11]  = 8.400e-1
  pft$clumping_factor[12:13] = 8.400e-1
  pft$clumping_factor[14:15] = clump_grass
  pft$clumping_factor[16]    = clump_grass
  pft$clumping_factor[17]    = clump_tree 
  
  
  onethird = 1./3.
  pft$leaf_turnover_rate = numeric(npft)
  pft$leaf_turnover_rate[1]          = 2.0
  pft$leaf_turnover_rate[2]          = 1.0
  pft$leaf_turnover_rate[3]          = 0.5
  pft$leaf_turnover_rate[4]          = onethird
  pft$leaf_turnover_rate[5]          = 2.0
  pft$leaf_turnover_rate[6]          = onethird
  pft$leaf_turnover_rate[7]          = onethird
  pft$leaf_turnover_rate[8]          = onethird
  pft$leaf_turnover_rate[9]          = 0.0
  pft$leaf_turnover_rate[10]         = 0.0
  pft$leaf_turnover_rate[11]         = 0.0
  pft$leaf_turnover_rate[12]         = 2.0
  pft$leaf_turnover_rate[13]         = 2.0
  pft$leaf_turnover_rate[14]         = 2.0
  pft$leaf_turnover_rate[15]         = 2.0
  pft$leaf_turnover_rate[16]         = 2.0
  pft$leaf_turnover_rate[17]         = 1.27
  
  C2B         = 2.
  
  sla_scale =  0.1 * C2B
  sla_inter =  2.4
  sla_slope = -0.46
  
  pft$SLA = numeric(npft)
  pft$SLA[1] = 22.7 
  pft$SLA[2] = 10.0**(sla_inter + sla_slope * log10(12.0/pft$leaf_turnover_rate[2] )) * sla_scale
  pft$SLA[3] = 10.0**(sla_inter + sla_slope * log10(12.0/pft$leaf_turnover_rate[3] )) * sla_scale
  pft$SLA[4] = 10.0**(sla_inter + sla_slope * log10(12.0/pft$leaf_turnover_rate[4] )) * sla_scale
  pft$SLA[5] = 22.0
  pft$SLA[6] =  6.0
  pft$SLA[7] =  9.0
  pft$SLA[8] = 10.0
  pft$SLA[9] = 30.0
  pft$SLA[10] = 24.2
  pft$SLA[11] = 60.0
  pft$SLA[12] = 22.0
  pft$SLA[13] = 22.0
  pft$SLA[14] = 22.7 # 10.0**(sla_inter + sla_slope * log10(12.0/pft$leaf_turnover_rate[14])) * sla_scale
  pft$SLA[15] = 22.7 # 10.0**(sla_inter + sla_slope * log10(12.0/pft$leaf_turnover_rate[15])) * sla_scale
  pft$SLA[16] = 22.7 
  pft$SLA[17] = 10.0**(sla_inter + sla_slope * log10(12.0/pft$leaf_turnover_rate[17])) * sla_scale
  
  
  pft$b1Bl_small = numeric(npft)
  pft$b1Bl_small[1:4]   = 0.0
  pft$b1Bl_small[5]     = 0.08
  pft$b1Bl_small[6]     = 0.024
  pft$b1Bl_small[7]     = 0.024
  pft$b1Bl_small[8]     = 0.0454
  pft$b1Bl_small[9]     = 0.0129
  pft$b1Bl_small[10]    = 0.048
  pft$b1Bl_small[11]    = 0.017
  pft$b1Bl_small[12:13] = 0.08
  pft$b1Bl_small[14:15] = 0.0
  pft$b1Bl_small[16]    = 0.0
  pft$b1Bl_small[17]    = 0.0
  
  pft$b2Bl_small = numeric(npft)
  pft$b2Bl_small[1:4]   = 0.0
  pft$b2Bl_small[5]     = 1.0
  pft$b2Bl_small[6]     = 1.899
  pft$b2Bl_small[7]     = 1.899
  pft$b2Bl_small[8]     = 1.6829
  pft$b2Bl_small[9]     = 1.7477
  pft$b2Bl_small[10]    = 1.455
  pft$b2Bl_small[11]    = 1.731
  pft$b2Bl_small[12:13] = 1.0
  pft$b2Bl_small[14:15] = 0.0
  pft$b2Bl_small[16]    = 0.0
  pft$b2Bl_small[17]    = 0.0
  
  pft$b1Bl_large     = pft$b1Bl_small
  pft$b2Bl_large     = pft$b2Bl_small
  
  pft$bleaf_adult = numeric(npft)
  
  pft$dbh_adult = numeric(npft)
  pft$dbh_adult [1:16] = 10.0
  pft$dbh_adult[17] = 1.81
  
  pft$min_dbh = numeric(npft)
  pft$min_dbh[1:16] = 0.1211817637
  pft$min_dbh[17] =  0.0478631184
  
  for(ipft in seq(1,npft)){
    if (pft$is_tropical[ipft]){
      if (allom < 2){
        pft$b1Bl_small [ipft] = exp(a1 + c1l * pft$b1Ht[ipft] + d1l * log(pft$rho[ipft]))
        aux = ( (a2l - a1) + pft$b1Ht[ipft] * (c2l - c1l) + log(pft$rho[ipft])* (d2l - d1l) ) * (1.0/log(dcrit))
        pft$b2Bl_small [ipft] = C2B * b2l + c2l * pft$b2Ht[ipft] + aux
        pft$b1Bl_large [ipft] = pft$b1Bl_small[ipft]
        pft$b2Bl_large [ipft] = pft$b2Bl_small[ipft]
        pft$bleaf_adult[ipft] = pft$b1Bl_large[ipft] / C2B * pft$dbh_adult[ipft] ** pft$b2Bl_large[ipft]
        
      } else if (allom == 2) {
        pft$b1Bl_small [ipft] = C2B * exp(nleaf[1]) * pft$rho[ipft] / nleaf[3]
        pft$b2Bl_small [ipft] = nleaf[2]
        pft$b1Bl_large [ipft] = pft$b1Bl_small[ipft]
        pft$b2Bl_large [ipft] = pft$b2Bl_small[ipft]
        pft$bleaf_adult[ipft] = pft$b1Bl_large[ipft] / C2B * pft$dbh_adult[ipft] ** pft$b2Bl_large[ipft]
        
      } else {
        pft$b1Bl_large [ipft] = 0.00873 * pft$SLA[3] / pft$SLA[ipft]
        pft$b2Bl_large [ipft] = 2.1360
        pft$bleaf_adult[ipft] = pft$b1Bl_large[ipft] / C2B * pft$dbh_adult[ipft] ** pft$b2Bl_large[ipft]
        pft$bleaf_sapling     = 0.02 * C2B * pft$SLA[3] / pft$SLA[ipft]
        pft$b2Bl_small [ipft] = log( pft$bleaf_adult[ipft] / pft$bleaf_sapling ) / log(pft$dbh_adult  [ipft] / pft$min_dbh[ipft] )
        pft$b1Bl_small [ipft] = pft$bleaf_adult[ipft] * C2B / pft$dbh_adult[ipft] ** pft$b2Bl_small[ipft]
        
      }
    } else {
      pft$bleaf_adult[ipft] = pft$b1Bl_large[ipft] / C2B * pft$dbh_adult[ipft] ** pft$b2Bl_large[ipft]
    }
  }
  
  pft$b1Bl_small[17] = pft$b1Bl_small[2]
  pft$b2Bl_small[17] = pft$b2Bl_small[2]
  pft$b1Bl_large[17] = 0.0643
  pft$b2Bl_large[17] = 1.6105
  
  allom_ed_param=pft
  return(allom_ed_param)
}
