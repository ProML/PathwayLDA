# Gibbs sampling iterations
# Parameters
alpha = .5
beta = prior_phi*.1 #asymmetric beta
num.iterations = 1000
burn_in = num.iterations/2
thinning_rate = 10
######################################################################
# Initialisation

cdt = matrix(0,nrow=num.docs,ncol=num.topics)
cvt = matrix(0,nrow=num.vocabs,ncol=num.topics)

z = vector(mode = "list", length = num.docs)
word.ids = vector(mode = "list", length = num.docs)
for (d in 1:num.docs) #each document
{
	z.d = vector()
	word.ids.d = vector()
	print(d)
	for (v in 1:num.vocabs) #each word
	{
	word.id = v
	freq = docs[d,v]
	w = 0
		while(w < freq)
		{
			w = w + 1
			topic.id = sample(num.topics,1,prob=prior_phi[v,]/sum(prior_phi[v,]))
			cdt[d,topic.id] = cdt[d,topic.id]+1
			#print(word.id)
			#print(prior_phi[word.id,topic.id])
			cvt[word.id,topic.id] = cvt[word.id,topic.id]+1 
			z.d = c(z.d,topic.id)
			word.ids.d = c(word.ids.d,word.id)
		}
	}
	z[[d]] = z.d
	word.ids[[d]] = word.ids.d
}

# Document-topic distributions
theta = prop.table(cdt+alpha,1)
	 
# Word-topic distributions
phi = prop.table(cvt+beta,2)




