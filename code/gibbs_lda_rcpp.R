library(Rcpp)
cppFunction('
using namespace Rcpp;
List gibbs_lda_rcpp(
double alpha, NumericMatrix beta, 
List z,List word_ids,NumericMatrix cdt,NumericMatrix cvt,NumericMatrix theta,NumericMatrix phi,
int num_iterations,int burn_in,int thinning_rate,
Function save_answers_r
)
{ 

int num_docs = cdt.nrow();
int num_topics = cdt.ncol();
int num_vocabs = cvt.nrow();
NumericVector ct_sum(num_docs);
NumericVector cv_sum(num_topics);
NumericVector beta_sum(num_topics);
for(int d=0;d<num_docs;d++)
{
	for(int t=0;t<num_topics;t++)
	{
		ct_sum(d) += cdt(d,t); 
	}
}
for(int t=0;t<num_topics;t++)
{
	for(int v=0;v<num_vocabs;v++)
	{
		cv_sum(t) += cvt(v,t); 
	}
}
for(int t=0;t<num_topics;t++)
{
	for(int v=0;v<num_vocabs;v++)
	{
		beta_sum(t) += beta(v,t); 
	}
}
for (int d=0;d<num_docs;d++) 
{ 
	NumericVector word_ids_d = word_ids(d);
	NumericVector topic_ids_d = z(d);
	for (int w=0;w<word_ids_d.size();w++)
	{
		word_ids_d(w)--;
		topic_ids_d(w)--;
	}
	word_ids(d) = word_ids_d;
	z(d) = topic_ids_d;
}
for (int iteration = 0; iteration < num_iterations;iteration++) 
{
   for (int d=0;d<num_docs;d++) 
   { 
Rcout << "Working on iteration:" << iteration << " document:" << d << "\\n";
		NumericVector word_ids_d = word_ids(d);
		NumericVector topic_ids_d = z(d);
		int num_words = word_ids_d.size();
     	for (int w=0;w<num_words;w++) 
	{
		//Rcout << "Working on iteration:" << iteration << " document:" << d << " word:" << w << "\\n";
		int word_id = word_ids_d(w);		
		int topic_old = topic_ids_d(w);
		
	// Decrement counts before computing equation (3)
        cdt(d,topic_old)--;
     	  cvt(word_id,topic_old)--; 
	  cv_sum(topic_old)--;
        ct_sum(d)--;

	// do multinomial sampling via cumulative method:
        NumericVector p(num_topics);
	  int topic_new;
	  int denom_theta = ct_sum(d) + num_topics * alpha;
        for (int t = 0; t<num_topics; t++) {
		phi(word_id,t) = (cvt(word_id,t) + beta(word_id,t))/(cv_sum(t)+beta_sum(t));
		theta(d,t) = (cdt(d,t) + alpha) / denom_theta; 
            p(t) = phi(word_id,t) * theta(d,t); //unnormalised
        }
        // cumulate multinomial parameters
        for (int t= 1; t<num_topics; t++) {
            p(t) += p(t - 1);
        }
        // scaled sample because of unnormalised p[]
        double u = R::runif(0,1) * p(num_topics - 1);
	  for (topic_new=0; topic_new<num_topics; topic_new++) {
            if (u < p(topic_new))
                break;
        }
       // add newly estimated z_i to count variables
        cvt(word_id,topic_new)++;
        cdt(d,topic_new)++;
        cv_sum(topic_new)++;
        ct_sum(d)++;
	  topic_ids_d(w) = topic_new;
		
		} // end w loop
	  z(d) = topic_ids_d;
	} // end d loop
		
	// Save samples
	save_answers_r(cdt,cvt,iteration);

} // end iteration loop
return Rcpp::List::create(Rcpp::Named("cvt_rcpp") = cvt,
                          Rcpp::Named("cdt_rcpp") = cdt,
			        Rcpp::Named("z_rcpp") = z
				  );
} //end function
')
# Gibbs sampling ends
save_answers <- function(cdt,cvt,iteration)
{
	write(t(cdt),file=paste("answers_cdt_",iteration,".txt",sep=""),sep = "\t",ncolumns=ncol(cdt))
	write(t(cvt),file=paste("answers_cvt_",iteration,".txt",sep=""),sep = "\t",ncolumns=ncol(cvt))
}
