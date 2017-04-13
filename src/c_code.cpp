#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List getNullHetParamC(List DataList, List GroupList, int ncluster) {
  int nsample = DataList.size();
  List d_list(nsample);
  List p_list(nsample);
  List rsid_list(nsample);
  List tmp1(1);
  List tmp2(1);
  for (int n=0; n<nsample; n++) {
    arma::mat x_null = DataList[n];
    int nr = x_null.n_rows;
    int nc = x_null.n_cols;
    arma::vec group_null = GroupList[n];
    // Ditance of random clusters
    arma::mat mu1 = arma::zeros(ncluster,nc);
    arma::mat mu2 = arma::zeros(ncluster,nc);
    arma::vec d = arma::zeros(ncluster,1);
    for(int i=0;i<ncluster;i++){
      arma::uvec id1 = find(group_null==(i+1));
      arma::vec rmean1 = arma::vectorise(arma::mean(x_null.rows(id1)));
      mu1.row(i) = rmean1.t();
      
      arma::uvec id2 = find(group_null!=(i+1));
      arma::vec rmean2 = arma::vectorise(arma::mean(x_null.rows(id2)));
      mu2.row(i) = rmean2.t();
      
      //Distance
      d(i) = sqrt(arma::sum(arma::pow(rmean1.t()-rmean2.t(),2)));
      // Rcout << "The value is " << d(i)<< std::endl;
    }
    d_list[n] = d;
    //Size of random clusters
    arma::vec p = arma::zeros(ncluster,1);
    for (int i=0;i<nr;i++){
      p(group_null(i)-1)++;
    }
    // arma::vec p_min = arma::min(p)*arma::ones(ncluster,1)/nr;
    arma::vec p_min = arma::min(p,nr-p)/nr;
    p_list[n] = p_min;
    // Sample index
    rsid_list[n] = (n+1)*arma::ones(ncluster,1);
  }
  return Rcpp::List::create(Rcpp::Named("d_list") = d_list,
                           Rcpp::Named("p_list") = p_list,
                           Rcpp::Named("rsid_list") = rsid_list);
}

// [[Rcpp::export]]
arma::vec cummaxC(arma::vec x){
  arma::vec y = arma::zeros(x.n_elem,1);
  for(int i=0;i<x.n_elem;i++){
    if(i==0){
      y(i) = x(i);
    }else{
      if(x(i)>y(i-1)){
        y(i) = x(i);
      }else{
        y(i) = y(i-1);
      }
    }
    
  }
  return y;
}

// [[Rcpp::export]]
arma::uvec descend2C(arma::vec x, arma::vec y){
  arma::uvec id;
  arma::uvec tmp_id1 = arma::stable_sort_index(y,"descend");
  arma::vec tmp_x = x(tmp_id1);
  arma::uvec tmp_id2 = arma::stable_sort_index(tmp_x,"descend");
  id = tmp_id1(tmp_id2);
  return id;
}

// [[Rcpp::export]]
List paretoFrontTest2C(arma::vec x, 
                       arma::vec y, 
                       arma::vec cluster_id,
                       arma::uvec cluster_size,
                       arma::vec x_null, 
                       arma::vec y_null,
                       arma::uvec rsid_null){
  int nsample = x_null.n_elem;
  arma::uvec index = arma::zeros<arma::uvec>(nsample);
  arma::vec flag = arma::zeros(nsample);
  arma::vec dom = arma::zeros(nsample);
  
  int level = 0;
  DataFrame df = DataFrame::create(Rcpp::Named("f1") = x_null,
                                   Rcpp::Named("f2") = y_null);
  
  while(sum(index == 0 && flag ==0)>0){
    level++;//index system starts from 1 (for R compatibility)
    arma::uvec id0 = arma::find(index == 0 && flag == 0);
    arma::vec f1 = x_null(id0);
    arma::vec f2 = y_null(id0);
    //     arma::uvec myorder = as<arma::uvec>(order(f1,
    //                                           f2,
    //                                           Rcpp::Named("decreasing") = true));
    // arma::uvec myorder = arma::sort_index(f1,"descend");
    arma::uvec myorder = descend2C(f1,f2);
    arma::vec f1_new = f1(myorder);
    arma::vec f2_new = f2(myorder);
    arma::vec cummax_f2 = cummaxC(f2_new);
    arma::vec cummax_f2_sub1 = cummax_f2.subvec(1,(cummax_f2.n_elem-1));
    arma::vec cummax_f2_sub2 = cummax_f2.subvec(0,(cummax_f2.n_elem-2));
    double eps = pow(10,-6);
    
    arma::vec cummax_f2_absdiff = abs(cummax_f2_sub1 - cummax_f2_sub2);
    arma::uvec is_not_duplicate = arma::zeros<arma::uvec>(cummax_f2_absdiff.n_elem+1);
    for(int i=0;i<is_not_duplicate.n_elem;i++){
      if(i==0){
        is_not_duplicate(i) = 1;
      }else{
        if(cummax_f2_absdiff(i-1)>eps){
          is_not_duplicate(i) = 1;
        }
      }
    }
    arma::uvec id1 = id0(myorder(arma::find(is_not_duplicate==1)));
    
    arma::uvec rsid = rsid_null(id1);
    
    // Only one sub-cluster from each sample is included in the null distribution
    arma::uvec uniq_rsid = arma::unique(rsid);
    for(int i = 0;i<uniq_rsid.n_elem;i++){
      arma::uvec tmp_id = arma::find(rsid_null==uniq_rsid(i));
      flag(tmp_id) = arma::ones<arma::vec>(tmp_id.n_elem);
      index(tmp_id) = arma::ones<arma::uvec>(tmp_id.n_elem)*level;
      arma::uvec tmp_id1 = arma::find(rsid==uniq_rsid(i));
      if(tmp_id1.n_elem==1){
        dom(id1(tmp_id1)) = arma::ones<arma::vec>(1);
      }else{
        arma::uvec id2 = id1(tmp_id1);
        dom(id2(find(id2==min(id2)))) = arma::ones<arma::vec>(1);
      }
    }
  }
  
  int nlevel = level;
  arma::vec pf_x_null = x_null(arma::find(dom==1));
  arma::vec pf_y_null = y_null(arma::find(dom==1));
  arma::uvec pf_rsid_null = rsid_null(arma::find(dom==1));
  arma::uvec pf_index = index(arma::find(dom==1));
  arma::vec pf_flag = flag(arma::find(dom==1));
  arma::vec pf_dom = dom(arma::find(dom==1));
  arma::vec pf_ref_pvals = arma::ones<arma::vec>(nlevel);
  
  // Assign p-values to Pareto sets
  for(int i=0;i<nlevel;i++){
    if(i==0){
      pf_ref_pvals(i) = 0;
    }else{
      // pf$ref.pvals[l] = length(which(pf$index<l))/length(pf$index)
      arma::uvec pf_l = arma::find(pf_index<(i+1));
      pf_ref_pvals(i) = pf_l.n_elem/double(pf_index.n_elem);
    }
  }
  
  // Significant test
  arma::mat pf_pval_obj_all = arma::zeros<arma::mat>(5,x.n_elem);
  double pval;
  for(int i=0;i<x.n_elem;i++){
    // which((pf$x.null<pf$x[i] & pf$y.null<=pf$y[i]) | (pf$x.null<=pf$x[i] & pf$y.null<pf$y[i]))
    arma::uvec tmp_sig_id = find((pf_x_null<x(i) && pf_y_null<=y(i))||(pf_x_null<=x(i) && pf_y_null<y(i)));
    if(tmp_sig_id.n_elem>0){
      int tmp_pf_id = arma::min(pf_index(tmp_sig_id));
      pval = pf_ref_pvals(tmp_pf_id-1);
    }else{
      pval = 1;
    }
    
    pf_pval_obj_all(0,i) = x(i);
    pf_pval_obj_all(1,i) = y(i);
    pf_pval_obj_all(2,i) = pval;
    pf_pval_obj_all(3,i) = cluster_size(i);
    if(cluster_size(i)<=(arma::sum(cluster_size)-cluster_size(i))){
      pf_pval_obj_all(4,i) = cluster_size(i);
    }else{
      pf_pval_obj_all(4,i) = arma::sum(cluster_size)-cluster_size(i);
    }
  }
  
  arma::uvec mdc_id = find(pf_pval_obj_all.row(2)==min(pf_pval_obj_all.row(2)));
  double pf_pval = pf_pval_obj_all(2,mdc_id(0));
  double pf_d_mdc = pf_pval_obj_all(0,mdc_id(0));
  double pf_p_mdc = pf_pval_obj_all(1,mdc_id(0));
  double pf_cluster_size_mdc = pf_pval_obj_all(3,mdc_id(0));
  double pf_cluster_eff_size_mdc = pf_pval_obj_all(4,mdc_id(0));
  double pf_cluster_id_mdc = cluster_id(mdc_id(0));
  
  
  
  return Rcpp::List::create(Rcpp::Named("x") = as<std::vector<double> >(wrap(x)),
                            Rcpp::Named("y") = as<std::vector<double> >(wrap(y)),
                            Rcpp::Named("cluster.id") = as<std::vector<double> >(wrap(cluster_id)),
                            Rcpp::Named("x.null") = as<std::vector<double> >(wrap(pf_x_null)),
                            Rcpp::Named("y.null") = as<std::vector<double> >(wrap(pf_y_null)),
                            Rcpp::Named("rsid") = as<std::vector<double> >(wrap(pf_rsid_null)),
                            Rcpp::Named("nlevel") = nlevel,
                            Rcpp::Named("index") = as<std::vector<double> >(wrap(pf_index)),
                            Rcpp::Named("flag") = as<std::vector<double> >(wrap(pf_flag)),
                            Rcpp::Named("dom") = as<std::vector<double> >(wrap(pf_dom)),
                            Rcpp::Named("ref.pvals") = as<std::vector<double> >(wrap(pf_ref_pvals)),
                            // Rcpp::Named("pval.obj.all") = as<std::vector<double> >(wrap(pf_pval_obj_all)),
                            Rcpp::Named("pval.obj.all") = wrap(pf_pval_obj_all),
                            Rcpp::Named("pval") = pf_pval,
                            Rcpp::Named("d.mdc") = pf_d_mdc,
                            Rcpp::Named("p.mdc") = pf_p_mdc,
                            Rcpp::Named("cluster.size.mdc") = pf_cluster_size_mdc,
                            Rcpp::Named("cluster.eff.size.mdc") = pf_cluster_eff_size_mdc,
                            Rcpp::Named("cluster.id.mdc") = pf_cluster_id_mdc);
}
// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R

*/
