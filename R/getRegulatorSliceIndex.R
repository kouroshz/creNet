getRegulatorSliceIndex = function(ents, rels, gene.ids)
{
   u.hyps = unique(rels$srcuid)
   child.uid = split(rels$trguid, rels$srcuid)
   child.sgn = numeric(nrow(rels))
   child.sgn[rels$type == "increase"] = 1
   child.sgn[rels$type == "decrease"] = -1
   child.sgn = split(child.sgn, rels$srcuid)
   ents.mRNA = ents[which(ents$type == "mRNA"), ]
   child.id = ents.mRNA$id[match(rels$trguid,ents.mRNA$uid)]
   data.slice = match(child.id, gene.ids)
   child.id = split(child.id,rels$srcuid)
   data.slice = split(data.slice,rels$srcuid)
   groups = rep(1:length(u.hyps), sapply(child.id, length))
    L = list(ents = ents, rels = rels, ents.mRNA = ents.mRNA, 
        u.hyps = u.hyps, child.uid = child.uid, child.sgn = child.sgn, 
        child.id = child.id, data.slice = data.slice, groups = groups)
    L
}

