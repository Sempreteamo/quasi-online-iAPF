normalizing_constant <- function(a0, P0, dt, ct, Tt, Zt, Ht, Gt, obs){
  fkf.obj <- fkf(a0, P0, dt, ct, Tt, Zt, Ht, Gt, yt = t(obs))
  fks.obj <- fks(fkf.obj)
  fkf.obj_Z <- fkf(a0, P0, dt, ct, Tt, Zt, Ht, Gt, yt = t(obs))$logLik
  
  return(fkf.obj_Z)
}

NC_estimate <- function(Z, fkf.obj_Z){
  print(exp(Z-fkf.obj_Z))
}
