
ico = 2
pftnum = PFT[ico]

b1Ca = param_all$b1Ca[pftnum]
b2Ca = param_all$b2Ca[pftnum]

b1Bl =  param_all$b1Bl_large[pftnum]
b2Bl =   param_all$b2Bl_large[pftnum]
n = nplant[ico]
dbh = DBH[ico]
sla = param_all$SLA[pftnum]
  
Bl = b1Bl/2 * (dbh**b2Bl)

CrownA = b1Ca*dbh**b2Ca
loclai = sla * Bl
dbh2ca = min (loclai, CrownA)

min(1.0, nplant * dbh2ca)
CA[ico]