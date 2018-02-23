getRegulatorSliceIndex = function(ents, rels, gene.ids)
{
   u.hyps = unique(rels$srcuid)
   ## split sorts the factor part which messes up the ordering.
   ## manually set the factor with desired levels
   ##child.uid = split(rels$trguid, rels$srcuid) 
   child.uid = split(rels$trguid, factor(rels$srcuid, levels = unique(rels$srcuid))) 
   child.sgn = numeric(nrow(rels))
   child.sgn[rels$type == "increase"] = 1
   child.sgn[rels$type == "decrease"] = -1
   ##same here
   ##child.sgn = split(child.sgn, rels$srcuid)
   child.sgn = split(child.sgn, factor(rels$srcuid, levels = unique(rels$srcuid)))
   ents.mRNA = ents[which(ents$type == "mRNA"), ]
   child.id = ents.mRNA$id[match(rels$trguid,ents.mRNA$uid)]
   data.slice = match(child.id, gene.ids)
   ## same here
   ##child.id = split(child.id,rels$srcuid)
   child.id = split(child.id,factor(rels$srcuid, levels = unique(rels$srcuid)))
   ##data.slice = split(data.slice,rels$srcuid)
   data.slice = split(data.slice,factor(rels$srcuid, levels = unique(rels$srcuid)))
   groups = rep(1:length(u.hyps), sapply(child.id, length))
    L = list(ents = ents, rels = rels, ents.mRNA = ents.mRNA, 
        u.hyps = u.hyps, child.uid = child.uid, child.sgn = child.sgn, 
        child.id = child.id, data.slice = data.slice, groups = groups)
    L
}

