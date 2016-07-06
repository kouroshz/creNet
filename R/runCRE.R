runCRE <- function(npp, npm, npz, nmp, nmm, nmz, nrp, nrm, nrz, nzp, nzm, nzz, method = c("Ternary","Quaternary","Enrichment","Fisher"))
{
  
  method = match.arg(method)
  
  if(method == 'Quaternary'){
    qR     = nrp + nrm + nrz
    qZero  = nzp + nzm + nzz
    nPlus  = npp + nmp + nrp + nzp
    nMinus = npm + nmm + nrm + nzm
    nZero  = npz + nmz + nrz + nzz
    
    ## Assume up-regulated
    qPlus  = npp + npm + npz
    qMinus = nmp + nmm + nmz
    score  = npp + nmm + nrp + nrm - (npm + nmp)
    pval.up   = QP_Pvalue(score = score, q_p = qPlus, q_m = qMinus, q_z = qZero,
                          q_r = qR, n_p = nPlus, n_m = nMinus, n_z = nZero)
    
    ## Assume down-regulated
    qPlus  = nmp + nmm + nmz
    qMinus = npp + npm + npz
    score  = nmp + npm + nrp + nrm - (npp + nmm)
    pval.down   = QP_Pvalue(score = score, q_p = qPlus, q_m = qMinus, q_z = qZero,
                            q_r = qR, n_p = nPlus, n_m = nMinus, n_z = nZero)
  } else if (method == 'Ternary'){
    qR     = 0
    qZero  = nzp + nzm + nzz
    nPlus  = npp + nmp + nzp
    nMinus = npm + nmm + nzm
    nZero  = npz + nmz + nzz
    
    ## Assume up-regulated
    qPlus  = npp + npm + npz
    qMinus = nmp + nmm + nmz
    score  = npp + nmm - (npm + nmp)
    pval.up   = QP_Pvalue(score = score, q_p = qPlus, q_m = qMinus, q_z = qZero,
                          q_r = qR, n_p = nPlus, n_m = nMinus, n_z = nZero)
    
    ## Assume down-regulated
    qPlus  = nmp + nmm + nmz
    qMinus = npp + npm + npz
    score  = nmp + npm - (npp + nmm)
    pval.down   = QP_Pvalue(score = score, q_p = qPlus, q_m = qMinus, q_z = qZero,
                            q_r = qR, n_p = nPlus, n_m = nMinus, n_z = nZero)
    
  } else if (method == 'Enrichment'){
    nrp    = npp + nmp + nrp
    nrm    = npm + nmm + nrm
    nrz    = npz + nmz + nrz
    
    qPlus  = 0
    qMinus = 0
    qR     = nrp + nrm + nrz
    qZero  = nzp + nzm + nzz
    
    nPlus  = nrp + nzp
    nMinus = nrm + nzm
    nZero  = nrz + nzz
    
    score  = nrp + nrm
    
    pval.up   = QP_Pvalue(score = score, q_p = qPlus, q_m = qMinus, q_z = qZero,
                          q_r = qR, n_p = nPlus, n_m = nMinus, n_z = nZero)
    pval.down = pval.up
  } else if (method == 'Fisher'){
    nrp    = npp + nmp + nrp
    nrm    = npm + nmm + nrm
    nrz    = npz + nmz + nrz
    
    M = matrix(0, nrow = 2, ncol = 2)
    M[1,1] = nrp + nrm
    M[1,2] = nrz
    M[2,1] = nzp + nzm
    M[2,2] = nzz
    
    pval.up = fisher.test(M, alternative = 'greater')$p.val
    pval.down = pval.up
    
  }
  
  pval = list(pval.up = pval.up, pval.down = pval.down)
  pval
}
