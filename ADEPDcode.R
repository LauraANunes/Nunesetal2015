
## This code is assuming investigator wants to downlist all species in the tree, one at a time.

library('phytools')
library('ape')
library('caper')

pext<-read.csv() #probabily of extinction of each species in the same order as tips of the tree!
phy<-read.tree('.tre') #get tree
tips<-phy$tip.label #extrat tip column 
tips<-data.frame(tips,c(1:length(tips))) #attribute index to locate species in the tree
nb.species<-length(tips) #length(tips) == number of species in the tree

colnames(tips)<-c('species','treeindex')

mat<-clade.matrix(phy)
brlen <- mat$edge.length #get branch lenghts in the same order as the p(ext) vector
nb.tips <- max(phy$edge) - phy$Nnode #number of tips
orphan <- nb.tips+1 ### indicate the basal branch (essentially nb.tips +1) this depends on whether the tree is rooted or unrooted.
nb.nodes <- phy$Nnode #number of nodes
nb.all <- nb.tips + nb.nodes #nodes+tips


#first find the parent nodes and corresponding children (sister pair) for whole tree

tab <- phy$edge #daughter of each edge 
tab <- data.frame(tab)  
names(tab)<-c('parent','child')
tab_child<-tab[with(tab,order(child)),] #order childs
tab_parent <- tab[with(tab,order(parent)),] #order parents
c1 <- rep(-1,nb.all) #child 1
c2 <- rep(-1,nb.all) #child 2
parent <- rep(-1,(nb.all-1)) #needs -1 to remove basalnode
for(i in 1:(nb.all-1))
{
  parent[i]<- c(tab_child[i,1])
}
parent <- append(parent,-1,after=(orphan-1))
c1[orphan] <- c(tab_parent[1,2]) #for child 1
### say p=(1,2) are the branches of the first sister pair, therefore p=3 is the branch of the parent
p<-3
for(m in 1:nb.all){
  c1[m+(nb.tips+1)] <- c(tab_parent[p,2])
  p <- p+2 #at every two steps (skip two sister pairs for next node)
}

c2[orphan] <-c(tab_parent[2,2]) #child 2
p <- 4
for(m in 1:nb.all){
  c2[m+(nb.tips+1)] <- c(tab_parent[p,2])
  p <- p+2 #at every two steps (skip two sister pairs for next node)
}

df<- cbind(c1,c2,parent)
df <- df[complete.cases(df),] #remove NAs

c1 <- as.vector(df[,1])
c2 <- as.vector(df[,2])
parent <- as.vector(df[,3])

# now obtain p(ext) of internal branches (parent branches) when species is downlisted
# add all p(ext)*branch lenghts values to get Expected Phylogenetic Diversity of whole tree


##current_epd == EDP with no downlisting 
parent_pext <- rep(NA,nb.nodes) ## at the moment is NA, will be recalculated next
all_pext <- c(pext,parent_pext) 

for(i in orphan:nb.all){
  sib1 <- c1[i]
  sib2 <- c2[i]
  parent_pext[i]<-pext[sib1]*pext[sib2] 
}
for(m in 1:60){ #repeats untill all all_pext have no NAs
  for(i in orphan:nb.all){
    sib1 <- c1[i]
    sib2 <- c2[i]
    all_pext[i]<-all_pext[sib1]*all_pext[sib2] 
  }
}

current_per <- brlen*(1-all_pext) # probability of survival of each branch
##make sure basal branch had p(ext) of 0.00

current_epd<- sum(as.numeric(current_per)) #current EXPECTED PHYLOGENETIC DIVERSITY (no downlisting)


level_ext<- c(0.00005,0.004,0.05,0.42,0.97,1) #probabilities of extinction according to IUCN50
runs <- c(1:nb.species) #1:number of species to downlist
adepd_epd <- rep(NA,nb.species)#number of species to downlist

for(j in 1:nb.species){
  cat(paste('',runs[j],'-m.  ',sep=''))
  ## use treeindex to locate where the species you want to downlist is located in the tree
  ## to downlist, you are moving the p(ext) one position below, there -1
  pext[j]<-level_ext[which(level_ext==pext[j])-1] #change the pext of species j to a lower level of extinction risk (level_ext)
  parent_pext <- rep(NA,nb.nodes)
  all_pext <- c(pext,parent_pext) #now recalculate all pext given change in pext[j]
  for(i in orphan:nb.all){
    sib1 <- c1[i]
    sib2 <- c2[i]
    parent_pext[i]<-pext[sib1]*pext[sib2] 
  }
  for(m in 1:60){ #repeats untill all all_pext have no NAs
    for(i in orphan:nb.all){
      sib1 <- c1[i]
      sib2 <- c2[i]
      all_pext[i]<-all_pext[sib1]*all_pext[sib2] 
    }
  }
  
  adepd.per <- brlen*(1-all_pext) # probability of survival of each branch
  ##make sure basal branch had p(ext) of 0.00
  
  adepd_epd[j]<- sum(as.numeric(adepd.per)) # EXPECTED PHYLOGENETIC DIVERSITY when species j is downlisted
  pext[j]<-level_ext[which(level_ext==pext[j])+1] #restore to original leveo of extinction risk of species 
}

ADEPD.score<-adepd_epd-current_epd #ADEPD=ADEPD_epd - current_epd as the gain in EPD due to the downlisting

res.frame <- data.frame(tips,adepd_epd,current_epd,ADEPD)  #species name, EPD of tree when species downlisted, EPD of tree with no downlisting, species ADEPD score

colnames(res.frame)<-c('species','ADEPD_EPD','CURRENT_EPD','ADEPD')  #frame with species index and ADEPD not divided by costs. 

#then divide by costs to get ADEPD/cost score
