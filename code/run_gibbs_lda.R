setwd("F:/Code/")
source('gibbs_lda_initialisation.R')
source('gibbs_lda_rcpp.R')

setwd("F:/Data/")
X = read.table("gene_expression_matrix_X.txt", header=FALSE, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)
K = read.table("gene_pathway_matrix_K.txt", header=FALSE, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)

docs = t(as.matrix(round(abs(X)*10)))
prior_phi = t(as.matrix(K))

num.topics = ncol(prior_phi)
num.docs = nrow(docs)
num.vocabs = ncol(docs)
num.words = as.vector(rowSums(docs))

setwd("F:/Results/")
result = gibbs_lda_rcpp(
alpha, beta, 
z,word.ids,cdt,cvt,theta,phi,
num.iterations,burn_in,thinning_rate,save_answers
)



