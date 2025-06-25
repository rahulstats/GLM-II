library(igraph)

Eng0 = initial_nu

Eng1 = matrix(NA, nrow=NROW(Eng0), ncol=2)

for (i in 1:nrow(Eng0)) {
  if(as.numeric(Eng0[i,3])<0) {Eng1[i,]= as.matrix(Eng0[i,c(2,1)])}
  else {Eng1[i,1:2]= as.matrix(Eng0[i,c(1,2)])}
}

## subsetting items
k1=8
Eng2 = subset(Eng1, Eng1[,1] %in% c(1:k1))
Eng2 = subset(Eng2, Eng2[,2] %in% c(1:k1))
## edges directions in vector
edges1= as.vector(t(apply(Eng2, 2, rev)))

g1=graph( edges=edges1, n=k1, directed=T )
pdf(file = "graph1.pdf")
plot(g1, vertex.size=25, edge.arrow.size=0.4, vertex.label.color="black",edge.color='black') #vertex.label.dist=1.5
dev.off()

KK=0
for (a in 1:k1) {
  KK[a] = sum(Eng2[,2]==a)
}








# read ascii file in r
dat <- read.csv("2023ANA.EVA", header = FALSE, skip = 1)
dat




mean(n_triad)

alpha= c(0.05,0.1,0.2,0.3,0.4,0.5)
size_mean= c(2.018,2.864,3.922,6.288,9.964,20.17)

plot(alpha, size_mean, type='l')



