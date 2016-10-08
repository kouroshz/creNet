### The following functions perform some statistics and post-processing:

########## Stat, pos-processing Start
## This function returns the network directly connected to a given node.
##
## Arguments:
##
## node     -    id of the node
##
## ents     -    entries of the network
##
## rels     -    relationships
##
## levels   -    boolean (include other level for proteins)
nodeNet = function(node, ents, rels, levels = F){
  node.ents = {}
  node.rels = {}
  
  ## Find the targets of the node
  targ = rels[which(rels[,2] == node),]
  node.ents = rbind(node.ents, ents[which(ents[,1] %in% targ[,3]), ])
  node.rels = rbind(node.rels, targ)
  
  ## Find the sources of the node
  src = rels[which(rels[,3] == node),]
  node.ents = rbind(node.ents, ents[which(ents[,1] %in% src[,2]), ])
  node.rels = rbind(node.rels, src)
  
  ## add the node
  node.ents = rbind(node.ents, ents[which(ents[,1] == node), ])
  
  ## consider other level
  if(levels){
    if(ents[which(ents[,1] == node), 4] == 'protein'){
      if(substr(node,nchar(node), nchar(node)) == '1'){
        node2 = paste(substr(node,1, (nchar(node) - 1)), '0', sep = '')
      }else{
        node2 = paste(substr(node,1, (nchar(node) - 1)), '1', sep = '')
      }
    }
    ## Find the targets of the node
    targ = rels[which(rels[,2] == node2),]
    node.ents = rbind(node.ents, ents[which(ents[,1] %in% targ[,3]), ])
    node.rels = rbind(node.rels, targ)
    
    ## Find the sources of the node
    src = rels[which(rels[,3] == node2),]
    node.ents = rbind(node.ents, ents[which(ents[,1] %in% src[,2]), ])
    node.rels = rbind(node.rels, src)
    
    ## add the node
    node.ents = rbind(node.ents, ents[which(ents[,1] == node2), ])
    
  }
  
  L = list(ents = unique(node.ents), rels = unique(node.rels))
  
}

## This function returns the network directly connected to a list of given node.
##
## Arguments:
##
## nodes    -    list of id of the nodes
##
## ents     -    entries of the network
##
## rels     -    relationships
##
## levels   -    boolean (include other level for proteins)
nodeNetList = function(nodes, ents, rels, levels = F){
  node.ents = {}
  node.rels = {}
  
  for(node in nodes){
    ## Find the targets of the node
    targ = rels[which(rels[,2] == node),]
    node.ents = rbind(node.ents, ents[which(ents[,1] %in% targ[,3]), ])
    node.rels = rbind(node.rels, targ)
    
    ## Find the sources of the node
    src = rels[which(rels[,3] == node),]
    node.ents = rbind(node.ents, ents[which(ents[,1] %in% src[,2]), ])
    node.rels = rbind(node.rels, src)
    
    ## add the node
    node.ents = rbind(node.ents, ents[which(ents[,1] == node), ])
    
    ## consider other level
    if(levels){
      if(ents[which(ents[,1] == node), 4] == 'protein'){
        if(substr(node,nchar(node), nchar(node)) == '1'){
          node2 = paste(substr(node,1, (nchar(node) - 1)), '0', sep = '')
        }else{
          node2 = paste(substr(node,1, (nchar(node) - 1)), '1', sep = '')
        }
      }
      ## Find the targets of the node
      targ = rels[which(rels[,2] == node2),]
      node.ents = rbind(node.ents, ents[which(ents[,1] %in% targ[,3]), ])
      node.rels = rbind(node.rels, targ)
      
      ## Find the sources of the node
      src = rels[which(rels[,3] == node2),]
      node.ents = rbind(node.ents, ents[which(ents[,1] %in% src[,2]), ])
      node.rels = rbind(node.rels, src)
      
      ## add the node
      node.ents = rbind(node.ents, ents[which(ents[,1] == node2), ])
      
    }
    
  }
  
  L = list(ents = unique(node.ents), rels = unique(node.rels))
  
}
