#include <Rcpp.h>
#include <vector>
#include "crn.cc"
using namespace Rcpp;
crn* grn=NULL;


// [[Rcpp::export]]
void rn_init(List input, int lsp, bool verbose = false){

  if(grn!=NULL) delete grn;
  
  grn=new crn(input,lsp,verbose);
  
  grn->gen_eqc();
  grn->gen_mgen();
}

// [[Rcpp::export]]
void rn_run_cs_gen(int n=20000000){
  
  if(grn==NULL){
    std::cout<<"Please init the reaction network (init_rn())"<<std::endl;
    return;
  }
  
  grn->cs_gen(n);
}

// [[Rcpp::export]]
void rn_run_c_cs_gen(){
  
  if(grn==NULL){
    std::cout<<"Please init the reaction network (init_rn())"<<std::endl;
    return;
  }
  
  grn->c_cs_gen();
}

// [[Rcpp::export]]
void rn_run_c_sso_gen(){
  
  if(grn==NULL){
    std::cout<<"Please init the reaction network (init_rn())"<<std::endl;
    return;
  }
  
  grn->c_sso_gen();
}


// [[Rcpp::export]]
NumericVector rn_run_sp2p(const NumericVector sp){
  
  if(grn==NULL){ 
    std::cout<<"Please init the reaction network (init_rn())"<<std::endl;
    return(0);
  }
  boost::dynamic_bitset<> temp =bitset_NV(sp,grn->lsp);
  boost::dynamic_bitset<> out = grn->sp2p(temp);
  return(bitget_NV(out));
}

// [[Rcpp::export]]
NumericVector rn_run_p2sp(const NumericVector p){
  
  if(grn==NULL){ 
    std::cout<<"Please init the reaction network (init_rn())"<<std::endl;
    return(0);
  }
  boost::dynamic_bitset<> temp =bitset_NV(p,grn->sp_p.size());
  boost::dynamic_bitset<> out = grn->p2sp(temp);
  return(bitget_NV(out));
}

// [[Rcpp::export]]
NumericVector rn_run_closure(const NumericVector p){
  
  if(grn==NULL){ 
    std::cout<<"Please init the reaction network (init_rn())"<<std::endl;
    return(0);
  }
  boost::dynamic_bitset<> temp =bitset_NV(p,grn->lsp);
  grn->closure(temp);
  return(bitget_NV(temp));
}

// [[Rcpp::export]]
bool rn_run_is_sep(const NumericVector p, bool conn =0){
  
  if(grn==NULL){
    std::cout<<"Please init the reaction network (init_rn())"<<std::endl;
    return(0);
  }
  boost::dynamic_bitset<> temp =bitset_NV(p,grn->lsp);
  return(grn->is_sep(temp));
}



// [[Rcpp::export]]
void rn_run_all_syn(bool no_bass =0){
  
  if(grn==NULL){ 
    std::cout<<"Please init the reaction network (init_rn())"<<std::endl;
    return;
  }
  
  grn->all_syn(no_bass); 
}

// [[Rcpp::export]]
void rn_run_cont_sep(bool conn=0){
  cout<<"entre"<<endl;
  if(grn==NULL){ 
    std::cout<<"Please init the reaction network (init_rn())"<<std::endl;
    return;
  }
  
  grn->cont_sep(conn);
}

// [[Rcpp::export]]
List rn_get(){
  
  if(grn==NULL){ 
    std::cout<<"Please init the reaction network (init_rn())"<<std::endl;
    return(0);
  }
  
  int n;
  
  n=grn->reac.size();
  std::vector<std::vector<double>  > reac(n) ;
  std::vector<std::vector<double>  > prod(n) ;
  for(int i=0;i<n;i++){
    reac[i]=bitget_v(grn->reac[i]);
    prod[i]=bitget_v(grn->prod[i]);
  }
  
  std::vector<std::vector<double>  > s_reac=grn->s_reac;
  std::vector<std::vector<double>  > s_prod=grn->s_prod;
  
  n=grn->sp_p.size();
  std::vector<std::vector<double>  > sp_p(n) ;
  for(int i=0;i<n;i++) sp_p[i]=bitget_v(grn->sp_p[i]);
  
  n=grn->sp_b.size();
  std::vector<std::vector<double>  > sp_b(n) ;
  for(int i=0;i<n;i++) sp_b[i]=bitget_v(grn->sp_b[i]);
  
  n=grn->r_p.size();
  std::vector<std::vector<double>  > r_p(n) ;
  for(int i=0;i<n;i++) r_p[i]=bitget_v(grn->r_p[i]);
  
  n=grn->r_b.size();
  std::vector<std::vector<double>  > r_b(n) ;
  for(int i=0;i<n;i++) r_b[i]=bitget_v(grn->r_b[i]);
  
  n=grn->rsp_b.size();
  std::vector<std::vector<double>  > rsp_b(n) ;
  for(int i=0;i<n;i++) rsp_b[i]=bitget_v(grn->rsp_b[i]);
  
  n=grn->p_b.size();
  std::vector<std::vector<double>  > p_b(n) ;
  for(int i=0;i<n;i++) p_b[i]=bitget_v(grn->p_b[i]);
  
  n=grn->conn.size();
  std::vector<std::vector<double>  > conn(n) ;
  for(int i=0;i<n;i++) conn[i]=bitget_v(grn->conn[i]);
  
  n=grn->mgen.size();
  std::vector<std::vector<std::vector<double> > > mgen(n);
  for(int i=0;i<n;i++){
    int m=grn->mgen[i].size();
    std::vector<std::vector<double> > mtmp(m);
    for(int j=0;j<m;j++) mtmp[j]=bitget_v(grn->mgen[i][j]);
    mgen[i]=mtmp;
  }
  
  n=grn->syn.size();
  std::vector<std::vector<double>  > syn(n);
  for(int i=0;i<n;i++) syn[i]=bitget_v(grn->syn[i]);
  
  n=grn->syn_clos.size();
  std::vector<std::vector<double>  > syn_clos(n);
  for(int i=0;i<n;i++) syn_clos[i]=bitget_v(grn->syn_clos[i]);
  
  n=grn->reac_part.size();
  std::vector<std::vector<double>  > reac_part(n) ;
  std::vector<std::vector<double>  > prod_part(n) ;
  for(int i=0;i<n;i++){
    reac_part[i]=bitget_v(grn->reac_part[i]);
    prod_part[i]=bitget_v(grn->prod_part[i]);
  }
  
  n=grn->clos_p.size();
  n=(n>1000000)? 1000000 : n;
  std::vector<std::vector<double>  > clos_sp(n);
  std::vector<std::vector<double>  > clos_p(n);
  for(int i=0;i<n;i++){ 
    clos_sp[i]=bitget_v(grn->clos_sp[i]);
    clos_p[i]=bitget_v(grn->clos_p[i]);
  }
  
  n=grn->c_clos_p.size();
  std::vector<std::vector<vector<std::vector<int> > > >c_clos_p(n);
  for(int i=0;i<n;i++){
    int m=grn->c_clos_p[i].size();
    c_clos_p[i]=vector<vector<vector <int> > >(m);
    int j=0;
    for (auto it=grn->c_clos_p[i].begin(); it!=grn->c_clos_p[i].end(); ++it){
      c_clos_p[i][j]=vector<vector <int> >(5);
      c_clos_p[i][j][0]=bitget_idx(it->first);
      c_clos_p[i][j][4]=bitget_idx(grn->p2sp(it->first));
      for(auto k=0;k<c_clos_p[i][j][0].size();k++) c_clos_p[i][j][0][k]++; 
      c_clos_p[i][j][1]=it->second[0];
      c_clos_p[i][j][2]=it->second[1];
      c_clos_p[i][j][3]=it->second[2];
      j++;
    }
  }
  
  n=grn->c_sso_p.size();
  std::vector<std::vector<vector<std::vector<int> > > >c_sso_p(n);
  for(int i=0;i<n;i++){
    int m=grn->c_sso_p[i].size();
    c_sso_p[i]=vector<vector<vector <int> > >(m);
    int j=0;
    for (auto it=grn->c_sso_p[i].begin(); it!=grn->c_sso_p[i].end(); ++it){
      c_sso_p[i][j]=vector<vector <int> >(5);
      c_sso_p[i][j][0]=bitget_idx(it->first);
      c_sso_p[i][j][4]=bitget_idx(grn->p2sp(it->first));
      for(auto k=0;k<c_sso_p[i][j][0].size();k++) c_sso_p[i][j][0][k]++; 
      c_sso_p[i][j][1]=it->second[0];
      c_sso_p[i][j][2]=it->second[1];
      c_sso_p[i][j][3]=it->second[2];
      j++;
    }
  }
  
  
  n=grn->ssmc_p.size();
  std::vector<std::vector<double>  > ssmc_sp(n);
  std::vector<std::vector<double>  > ssmc_p(n);
  for(int i=0;i<n;i++){
    ssmc_sp[i]=bitget_v(grn->ssmc_sp[i]);
    ssmc_p[i]=bitget_v(grn->ssmc_p[i]);
  }
  
  n=grn->ssmo_p.size();
  std::vector<std::vector<double>  > ssmo_sp(n);
  std::vector<std::vector<double>  > ssmo_p(n);
  for(int i=0;i<n;i++){
    ssmo_sp[i]=bitget_v(grn->ssmo_sp[i]);
    ssmo_p[i]=bitget_v(grn->ssmo_p[i]);
  }

  List L=List::create();
  L.push_back(reac,"reac");
  L.push_back(prod,"prod");
  L.push_back(s_reac,"s_reac");
  L.push_back(s_prod,"s_prod");
  L.push_back(sp_p,"sp_p");
  L.push_back(sp_b,"sp_b");
  L.push_back(r_p,"r_p");
  L.push_back(r_b,"r_b");
  L.push_back(rsp_b,"rsp_b");
  L.push_back(p_b,"p_b");
  L.push_back(IntegerVector(std::begin(grn->x_r_p),std::end(grn->x_r_p)),"x_r_p");
  L.push_back(conn,"conn");
  L.push_back(mgen,"mgen");
  L.push_back(syn,"syn");
  L.push_back(syn_clos,"syn_clos");
  L.push_back(reac_part,"reac_part");
  L.push_back(prod_part,"prod_part");
  L.push_back(IntegerVector(std::begin(grn->xsyn_eqc),std::end(grn->xsyn_eqc)),"xsyn_eqc");
  L.push_back(clos_sp,"clos_sp");
  L.push_back(clos_p,"clos_p");
  L.push_back(c_clos_p,"c_clos_p");
  L.push_back(c_sso_p,"c_sso_p");
  L.push_back(grn->cnt_clos,"cnt_clos");
  L.push_back(IntegerVector(std::begin(grn->n_clos_sp),std::end(grn->n_clos_sp)),"n_clos_sp");
  L.push_back(IntegerVector(std::begin(grn->n_clos_p),std::end(grn->n_clos_p)),"n_clos_p");
  L.push_back(IntegerVector(std::begin(grn->n_cont_clos_p),std::end(grn->n_cont_clos_p)),"n_cont_clos_p");
  L.push_back(IntegerVector(std::begin(grn->n_sep_clos_p),std::end(grn->n_sep_clos_p)),"n_sep_clos_p");
  L.push_back(IntegerVector(std::begin(grn->n_uncon_clos_p),std::end(grn->n_uncon_clos_p)),"n_uncon_clos_p");
  L.push_back(bitget_v(grn->rsm),"rsm");
  L.push_back(bitget_v(grn->ssm),"ssm");
  L.push_back(ssmc_sp,"ssmc_sp");
  L.push_back(ssmc_p,"ssmc_p");
  L.push_back(ssmo_sp,"ssmo_sp");
  L.push_back(ssmo_p,"ssmo_p");
  
  return(L);
}
/*** R

*/
