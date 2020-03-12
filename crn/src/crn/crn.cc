#include "crn.h"

using namespace Rcpp;
using namespace std;
using namespace boost;


// void reorder(vector<dynamic_bitset<> >& vA, vector<size_t>& vOrder)  
// {   
//   assert(vA.size() == vOrder.size());
//   
//   // for all elements to put in place
//   for( int i = 0; i < vA.size() - 1; ++i )
//   { 
//     // while the element i is not yet in place 
//     while( i != vOrder[i] )
//     {
//       // swap it with the element at its final place
//       int alt = vOrder[i];
//       swap( vA[i], vA[alt] );
//       if(i==0) cout<<vA[0]<<" reordered"<<endl;
//       swap( vOrder[i], vOrder[alt] );
//     }
//   }
// }

template< class T >
void reorder(
    std::vector<T> & unordered, 
    std::vector<size_t> const & index_map, 
    std::vector<T> & ordered)
{
  // copy for the reorder according to index_map, because unsorted may also be
  // sorted
  std::vector<T> copy = unordered;
  ordered.resize(index_map.size());
  for(int i = 0; i<index_map.size();i++)
  {
    ordered[i] = copy[index_map[i]];
  }
}

vector<size_t> sort_indexes(const vector<dynamic_bitset <> > &v) {
  
  // initialize original index locations
  vector<size_t> idx(v.size());
  iota(idx.begin(), idx.end(), 0);
  
  // sort indexes based on comparing values in v
  sort(idx.begin(), idx.end(),[&v](size_t i1, size_t i2) {return v[i1].count() < v[i2].count();});
  
  return idx;
}

bool comp_reac(dynamic_bitset<> r1,dynamic_bitset<> r2){ 
  return(r1.count() < r2.count());
}


void v2bitset(dynamic_bitset<> & b, const vector<double> &idx){
  for(int i=0;i<idx.size();i++) b.set(idx[i]-1,1);
}


void vset(vector<double> &b, const vector<double> &idx, const vector<double> &s){
  for(int i=0;i<idx.size();i++) b[idx[i]-1]=s[i];
}

vector<double> bitget_v(const dynamic_bitset<> & b, int offset = 1){
  
  vector<double> v;

  dynamic_bitset<>::size_type i = b.find_first();
  
  while(i!=dynamic_bitset<>::npos){
    v.push_back(i+offset);
    i=b.find_next(i);
  }
  
  return(v);
}

vector<int> bitget_idx(const dynamic_bitset<> & b){
  
  vector<int> v;
  
  dynamic_bitset<>::size_type i = b.find_first();
  
  while(i!=dynamic_bitset<>::npos){
    v.push_back(i);
    i=b.find_next(i);
  }
  
  return(v);
}

NumericVector bitget_NV(dynamic_bitset<> & b){
  
  NumericVector v;
  
  dynamic_bitset<>::size_type i = b.find_first();
  
  while(i!=dynamic_bitset<>::npos){
    v.push_back(i+1);
    i=b.find_next(i);
  }
  
  
  return(v);
}


dynamic_bitset <> bitset_NV(const NumericVector & v, int tp){
  
  dynamic_bitset <> r(tp);
  for(int i =0; i<v.size();i++){
    if(v[i]<1 || v[i] >tp){ 
      cout<<"index out of set range"<<endl;
      break;
    }
    else r.set(v[i]-1,1);
  }
  return(r);
}

 
// int bs_comp(dynamic_bitset <> s1, dynamic_bitset <> s2){
// 
//   int c1=s1.count(), c2=s2.count();
//   if(c1<c2) return(-1);
//   else if(c1>c2) return(1);
//   else if(s1==s2) return(0);
//   else if(s1<s2) return(-1);
//   else return(1);
// 
//   }

//Which a greater than b
dynamic_bitset <> which_gt(const vector<double> a,const vector<double> b){
  dynamic_bitset <> r(a.size());
  for(int i=0;i<a.size();i++) if(a[i]>b[i]) r.set(i);
  return(r);
}



//Empty Constructor
crn::crn(int sp, int rc){
  
  lsp=sp;
  lrn=rc;
  reac=vector<dynamic_bitset<> >(lrn);
  prod=vector<dynamic_bitset<> >(lrn);
  for(int i=0; i<lrn; i++){
    reac.at(i)=dynamic_bitset<>(lsp);
    prod.at(i)=dynamic_bitset<>(lsp);
  }  
}

//Constructor with List argument typ crn
crn::crn(List input, int sp, bool verbose){
  
  lsp=sp;
  vector<string> names = as<List>(input[0]).names();
  lrn=input.size();
  reac=vector<dynamic_bitset<> >(lrn);
  prod=vector<dynamic_bitset<> >(lrn);
  s_reac=vector<vector<double> >(lrn);
  s_prod=vector<vector<double> >(lrn);

  
  for(int i=0; i<lrn; i++){
    reac.at(i)=dynamic_bitset<>(lsp);
    prod.at(i)=dynamic_bitset<>(lsp);
    s_reac[i]=vector<double> (lsp);
    s_prod[i]=vector<double> (lsp);

    
  }
  
  
  for(int i = 0; i<input.size();i++){
    List rni = input[i];
    
    v2bitset(reac[i],rni[0]);
    v2bitset(prod[i],rni[1]);
    
    vset(s_reac[i],rni[0],rni[2]);
    vset(s_prod[i],rni[1],rni[3]);

    if(verbose) for(int j = 0; j<rni.size(); j++){
      cout<<names[j]<<":";
      vector<double> st = rni[j];
      for(int k=0; k<st.size(); k++){
        cout<<" "<<st.at(k);
      }
      cout<<endl; 
    } 
  }
}

void crn::closure(dynamic_bitset<> & sp){

  vector <int> n_reac(lrn);
  
  for(int i=0;i<lrn;i++){
    n_reac.at(i)=i;
  }
  
  int i=0;
  bool flag=0;
  while(n_reac.size()>0){
    
    if(reac[n_reac[i]].is_subset_of(sp)){
      sp=sp|prod[n_reac[i]];
      n_reac.erase(n_reac.begin()+i);
      flag=1;
    }else i++;
    if(i>=n_reac.size()){
      if(flag){ 
        i=0; 
        flag=0;
      }else break;
    }
  }
}


void crn::clos_one(const dynamic_bitset<> & sp, dynamic_bitset<> & res){
  
  res=sp;
  for(int i=0;i<lrn;i++) if(reac[i].is_subset_of(sp)) res=res|prod[i];
    
}

bool crn::is_rsm(const dynamic_bitset<> &sp){
  dynamic_bitset<> r(lsp);
  dynamic_bitset<> p(lsp);
  for(int i=0;i<lrn;i++){ 
    if(reac[i].is_subset_of(sp)){ 
      r|=reac[i];
      p|=prod[i];
    }
  }
  return(r.is_subset_of(p));
}

bool crn::is_ssm(const dynamic_bitset<> &sp){
  
  
  dynamic_bitset<> r(lsp);
  dynamic_bitset<> p(lsp);
  for(int i=0;i<lrn;i++){ 
    if(reac[i].is_subset_of(sp)){
      for(int j=0;j<lsp;j++){
        double d=s_prod[i][j]-s_reac[i][j];
        if(d<0) r.set(j,1);
        else if(d>0) p.set(j,1);
      }
    }
  }
  return(r.is_subset_of(p));
}

void crn::gen_eqc(){
  
  vector<dynamic_bitset<> > c_reac = reac;
  
  for(int i=0;i<lrn;i++) closure(c_reac[i]);
  
  int xeqc=0;
  x_r_p=vector<int>(lrn);
  x_r_p.assign(lrn,-1);
  
  for(int i=0;i<lrn;i++){
    if(x_r_p[i]>=0) continue;
    x_r_p[i]=xeqc;
    sp_b.push_back(c_reac[i]);
    
    for(int j=i+1;j<lrn;j++){
      
      if(c_reac[i]==c_reac[j]) x_r_p[j]=xeqc;
    }
    xeqc++;
  }
  
  
  sp_p = vector<dynamic_bitset<> >(xeqc);
  for(int i=0;i<xeqc;i++) sp_p[i]=dynamic_bitset<> (lsp);
  for(int i=0;i<lrn;i++) sp_p[x_r_p[i]]|=reac[i]|prod[i];
  
  r_p = vector<dynamic_bitset<> >(xeqc);
  for(int i=0;i<xeqc;i++) r_p[i]=dynamic_bitset<> (lrn);
  for(int i=0;i<lrn;i++) r_p[x_r_p[i]][i]=1;
  
  p_b=vector<dynamic_bitset<> >(sp_b.size());
  r_b=vector<dynamic_bitset<> >(sp_b.size());
  rsp_b=vector<dynamic_bitset<> >(sp_b.size());
  sp_sp_b=vector<dynamic_bitset<> >(sp_b.size());
  sn_sp_b=vector<dynamic_bitset<> >(sp_b.size());
  for(int i=0;i<p_b.size();i++){ 
    p_b[i]=sp2p(sp_b[i]);
    r_b[i]=dynamic_bitset<> (lrn);
    rsp_b[i]=dynamic_bitset<> (lsp);
    sp_sp_b[i]=dynamic_bitset<> (lsp);
    sn_sp_b[i]=dynamic_bitset<> (lsp);
    
    for(int j=p_b[i].find_first();j!=dynamic_bitset<>::npos;j=p_b[i].find_next(j)){ 
      r_b[i]|=r_p[j];
      for(int k=r_p[j].find_first();k!=dynamic_bitset<>::npos;k=r_p[j].find_next(k)){ 
        rsp_b[i]|=reac[k];
        sp_sp_b[i]|=which_gt(s_reac[k],s_prod[k]);
        sn_sp_b[i]|=which_gt(s_prod[k],s_reac[k]);
      }
    }
  }
  conn=vector<dynamic_bitset<> >(sp_b.size());
  
  for(int i=0;i<conn.size();i++) conn[i]=dynamic_bitset<> (conn.size());
  
  for(int i=0;i<conn.size();i++){
    for(int j=0;j<conn.size();j++){
      // connectivity conditions 
      //if(j!=i && p_b[i][j]==0 && rsp_b[j].intersects(sp_b[i])){  
      //if(j!=i && p_b[i][j]==0 && sp_b[j].intersects(sp_b[i])){   
      //if(j!=i && p_b[i][j]==0){   
      //Hasse condition (all connected)
      if(true){
        conn[i].set(j,1); 
        conn[j].set(i,1);
      }
    }
  }
}

dynamic_bitset<> crn::conn_b(const dynamic_bitset <> &s){
  dynamic_bitset<> c(s.size());
  for(int i=s.find_first();i!=dynamic_bitset<>::npos;i=s.find_next(i)) c|=conn[i];
  c-=s;
  return(c);
}

dynamic_bitset<> crn::contrib_b(const dynamic_bitset <> &s){
  dynamic_bitset<> p(lsp);
  dynamic_bitset<> n(lsp);
  for(int i=s.find_first();i!=dynamic_bitset<>::npos;i=s.find_next(i)){
    p|=sp_sp_b[i];
    n|=sn_sp_b[i];
  } 
  n-=p;
  
  dynamic_bitset<> c=s;
  if(n.empty()) return(c.reset());
  c.flip();
  dynamic_bitset<> pp(lsp);
  for(int i=c.find_first();i!=dynamic_bitset<>::npos;i=c.find_next(i)){
    if(!sp_sp_b[i].intersects(n)) c.set(i,0);
    pp|=sp_sp_b[i];
  }
  
  if((n-pp).empty()) c.reset();
  return(c);
}


bool crn::are_conn(const dynamic_bitset<> s1 , const dynamic_bitset<> s2){ return(s1.intersects(s2));}


void crn::gen_mgen(){
  
  
  for(int i=0;i<r_p.size();i++){
    
    vector<dynamic_bitset<> > v(r_p[i].count());
    int k=0;
    for(int j=r_p[i].find_first();j!=dynamic_bitset<>::npos;j=r_p[i].find_next(j)) v[k++]=reac[j];
    
    stable_sort(v.begin(), v.end(), comp_reac);
    
    
    for(int j=0;j< v.size()-1; j++){
      
      for(int k=j+1;k< v.size();){
        
        if(v[j].is_subset_of(v[k])) v.erase(v.begin()+k);
        else k++;
      }
      
    }
    mgen.push_back(v);
  }
  
}


void crn::cs_gen(int N=20000000){
  
  
  N_clos=N;
  clos_sp=vector<dynamic_bitset<> >(N);
  clos_p=vector<dynamic_bitset<> >(N);
  
  cnt_clos=0;
  
  dynamic_bitset<> pb(r_p.size());
  pb.flip();
  
  for(int i=0; i<r_p.size(); i++){
    if(!r_cs_gen(pb,sp_b[i],i)) break;
    pb.set(i,0);
  }
  
  clos_sp.erase(clos_sp.begin()+cnt_clos,clos_sp.end());
  clos_p.erase(clos_p.begin()+cnt_clos,clos_p.end());
  
  n_clos_sp=vector<int>(clos_sp.size());
  n_clos_p=vector<int>(clos_sp.size());
  
  vector <size_t> ind = sort_indexes(clos_p);
  
  reorder(clos_p,ind,clos_p);
  // reorder(clos_sp,ind);
  reorder(clos_sp,ind,clos_sp);
  
  
  for(int i=0;i<clos_sp.size();i++){
    n_clos_sp[i]=clos_sp[i].count();
    n_clos_p[i]=clos_p[i].count();
  }
}

bool crn::r_cs_gen(dynamic_bitset<> pb,dynamic_bitset<> sp, int o){
// pb: cadidates to be added
// sp: added species
  dynamic_bitset<> csp=sp;
  closure(csp);
  dynamic_bitset<> ib = sp2p(csp); 
  
  if(!ib.is_subset_of(pb)) return(1); //contained bassic
  clos_sp[cnt_clos]=csp;
  clos_p[cnt_clos]=sp2p(csp);
  cnt_clos++;
  if(cnt_clos>=N_clos) return(0); 
  o++;
  
  for(int i=o; i<r_p.size(); i++){
    if(ib[i]==1) continue;
    if(p_b[i].is_subset_of(pb)) if(!r_cs_gen(pb,csp|sp_b[i],i)) return(0);
    pb.set(i,0);
  }
return(1);
}

void crn::c_cs_gen() {
  int nb=r_b.size();
  c_clos_p = vector <map <dynamic_bitset <>, vector<vector<int> > > >(nb);
  dynamic_bitset<> tmp_sp;
  for(int i=0;i<nb;i++) c_clos_p[i] = map <dynamic_bitset <>, vector<vector<int> > >();
  for(int i=0;i<nb;i++){ 
    c_clos_p[p_b[i].count()-1][p_b[i]]=vector<vector<int> >(3);
    if(is_ssm(sp_b[i])){
      ssmc_sp.push_back(sp_b[i]);
      ssmc_p.push_back(p_b[i]);
    }
  }
  
  for(int l=0;l<nb;l++){
    cout<<"level="<<l<<endl;
    int k_i=1;
    for (auto k=c_clos_p[l].cbegin(); k != c_clos_p[l].cend(); ++k){
      dynamic_bitset<> c=conn_b((*k).first);
      
      // dynamic_bitset<> c=(*k).first;
      // c.flip();
      for(int j=c.find_first();j!=dynamic_bitset<>::npos;j=c.find_next(j)){
        dynamic_bitset<> tmp=(*k).first;
        tmp.set(j,1);
        tmp_sp=p2sp(tmp);
        closure(tmp_sp);
        tmp=sp2p(tmp_sp);
        
        map<dynamic_bitset <>, vector<vector<int> > >::iterator it=c_clos_p[tmp.count()-1].find(tmp);
        if(it==c_clos_p[tmp.count()-1].end()){ 
          c_clos_p[tmp.count()-1][tmp]={{l+1},{k_i},{j+1}};
          if(is_ssm(tmp_sp)){
            ssmc_sp.push_back(tmp_sp);
            ssmc_p.push_back(tmp);
          }
        }
        else{
          (*it).second[0].push_back(l+1);
          (*it).second[1].push_back(k_i);
          (*it).second[2].push_back(j+1);
        }
      }
       k_i++;
    }
  }
}

void crn::c_sso_gen() {
  int nb=r_b.size();
  c_sso_p = vector <map <dynamic_bitset <>, vector<vector<int> > > >(nb);
  dynamic_bitset<> tmp_sp;
  for(int i=0;i<nb;i++) c_sso_p[i] = map <dynamic_bitset <>, vector<vector<int> > >();
  for(int i=0;i<nb;i++){ 
    c_sso_p[p_b[i].count()-1][p_b[i]]=vector<vector<int> >(3);
  }
  
  for(int l=0;l<nb;l++){
    cout<<"level="<<l<<endl;
    int k_i=1;
    for (auto k=c_sso_p[l].cbegin(); k != c_sso_p[l].cend(); ++k){
      tmp_sp=p2sp((*k).first);
      if(is_ssm(tmp_sp)){
        
        ssmo_sp.push_back(tmp_sp);
        ssmo_p.push_back((*k).first);
        
        dynamic_bitset<> c=conn_b((*k).first);
        
        for(int j=c.find_first();j!=dynamic_bitset<>::npos;j=c.find_next(j)){
          dynamic_bitset<> tmp=(*k).first;
          tmp.set(j,1);
          tmp_sp=p2sp(tmp);
          closure(tmp_sp);
          tmp=sp2p(tmp_sp);
          
          map<dynamic_bitset <>, vector<vector<int> > >::iterator it=c_sso_p[tmp.count()-1].find(tmp);
          if(it==c_sso_p[tmp.count()-1].end()){ 
            c_sso_p[tmp.count()-1][tmp]={{l+1},{k_i},{j+1}};
          }
          else{
            (*it).second[0].push_back(l+1);
            (*it).second[1].push_back(k_i);
            (*it).second[2].push_back(j+1);
          }
        }
      }
      else{
        
        dynamic_bitset<> c=contrib_b((*k).first);
        
        for(int j=c.find_first();j!=dynamic_bitset<>::npos;j=c.find_next(j)){
          dynamic_bitset<> tmp=(*k).first;
          tmp.set(j,1);
          tmp_sp=p2sp(tmp);
          closure(tmp_sp);
          tmp=sp2p(tmp_sp);
          
          map<dynamic_bitset <>, vector<vector<int> > >::iterator it=c_sso_p[tmp.count()-1].find(tmp);
          if(it==c_sso_p[tmp.count()-1].end()){ 
            c_sso_p[tmp.count()-1][tmp]={{l+1},{k_i},{j+1}};
          }
          else{
            (*it).second[0].push_back(l+1);
            (*it).second[1].push_back(k_i);
            (*it).second[2].push_back(j+1);
          }
        }

      }
      k_i++;
    }
  }
}



//--------------------------

void crn::syn_gen(dynamic_bitset<> sp){
  
  vector<int> xs;
  vector<int> xp;
  
  for(int i=sp.find_first();i!=dynamic_bitset<>::npos;i=sp.find_next(i)) xs.push_back(i);
  
  vector<dynamic_bitset<> > sspp;
  
  dynamic_bitset<> b(xs.size());
  
  for(int i=0; i< sp_p.size();i++){
    
    if(sp.intersects(sp_p[i])){
      xp.push_back(i);
      for(int j=0;j<xs.size();j++) b[j]=sp_p[i][xs[j]];
      sspp.push_back(b);
      b.reset();
    }
  }
  
  dynamic_bitset<> p(xp.size());
  for(int i=0; i<xp.size(); i++) r_syn_gen(p,i,sspp,xp);
  
  
}

//---------------------------

void crn::r_syn_gen( dynamic_bitset<> &p, 
                     int o,
                     vector<dynamic_bitset<> > & sspp,
                     vector<int> &xp){
  
  //bool flag = xp[o]==9 && p.count()==0;
  dynamic_bitset<> u(sspp[0].size());
  p[o]=1;
  //if(flag) cout<<"r_syn_gen, p="<<p<<endl;
  if(p.count()>1) for(int i=p.find_first();i!=dynamic_bitset<>::npos;i=p.find_next(i)){
    p[i]=0;
    u.reset();
    for(int j=p.find_first();j!=dynamic_bitset<>::npos;j=p.find_next(j)) u|=sspp[j];
    p[i]=1;
    if(sspp[i].is_subset_of(u)){ 
      p[o]=0;
      return;
    }
    
  }
  
  
  u|=sspp[o];
  //if(flag) cout<<"u="<<u<<" all(u)="<<(int) u.all()<<endl;
  if(u.all()){
    
    dynamic_bitset<> c_syn(sp_p.size());
    
    for(int j=p.find_first();j!=dynamic_bitset<>::npos;j=p.find_next(j)) c_syn[xp[j]]=1;
    //if(flag) cout<<"c_syn="<<c_syn<<endl;
    syn.push_back(c_syn);
    p[o]=0;
    return;
  }
  
  for(int i=o+1; i<xp.size(); i++) r_syn_gen(p,i,sspp,xp);
  p[o]=0;
}


//---------------------------

void crn::all_syn(bool no_bass=1){
  
  for(int i=0;i<mgen.size();i++){ 
    
    for(int j=0;j<mgen[i].size();j++) {
      cout<<"part="<<i+1<<" n_gen="<<j+1<<" gen size="<<mgen[i][j].count()<<endl;
      if(mgen[i][j].count()>0) syn_gen(mgen[i][j]);
    }
  }
  for(int j=0;j< syn.size()-1; j++){
    for(int k=j+1;k< syn.size();){
      if(syn[j]==syn[k]) syn.erase(syn.begin()+k);
      else k++;
    }
  }
  
  cout<<syn.size()<<" synergies"<<endl;
  cout<<sp_b.size()<<" bassics"<<endl;
  vector<dynamic_bitset<> > sp_syn(syn.size());
  
  //species of the proto-synergies
  for(int i=0;i<syn.size();i++){
    int j=syn[i].find_first();
    sp_syn[i]=sp_p[j];
    for(int k=syn[i].find_next(j);k!=dynamic_bitset<>::npos;k=syn[i].find_next(k)) sp_syn[i]|=sp_p[k];
   }
  
  //partiton reactions
  for(int i=0; i<syn.size();i++){
    dynamic_bitset<> syn_p_sp;
    clos_one(sp_syn[i],syn_p_sp);
    dynamic_bitset<> p=sp2p(syn_p_sp);
    p-=syn[i];
    if(p.any()){
      reac_part.push_back(syn[i]);
      prod_part.push_back(p);
    }
  }
  
  //closure (species) of the proto-synergies
  for(int i=0;i<syn.size();i++) closure(sp_syn[i]);
  
  
  if(no_bass) for(int i=0; ;){
    next_case:
    if(i>=sp_syn.size()) break;
     for(int j=0;j<sp_b.size();j++) if(sp_b[j]==sp_syn[i]){ 
      sp_syn.erase(sp_syn.begin()+i);
      syn.erase(syn.begin()+i);
       goto next_case;
    }
    i++;  
  }
  
  
  //equivalence classes of synergies
  int xeqc=0;
  int xnsp=0;
  xsyn_eqc=vector<int>(syn.size());
  xsyn_eqc.assign(syn.size(),-1);

  //cout<<"486="<<sp_syn[485]<<endl;
  //cout<<"512="<<sp_syn[511]<<endl;


  for(int i=0;i<syn.size();i++){
    //if(i==511) cout<<"tag with="<<xsyn_eqc[i]<<endl;
    if(xsyn_eqc[i]>=0) continue;
    xsyn_eqc[i]=xeqc;
    syn_clos.push_back(sp_syn[i]);

    if(!is_sep(sp_syn[i])) xnsp++;

    for(int j=i+1;j<syn.size();j++){
      if(xsyn_eqc[j]==-1 && sp_syn[i]==sp_syn[j]) xsyn_eqc[j]=xeqc;
      //if(j==485 && xsyn_eqc[j]>0) cout<<"tag with="<<xeqc<<endl;
      //if(j==511 && xsyn_eqc[j]>0) cout<<"tag with="<<xeqc<<endl;
    }
    xeqc++;
  }

  
  cout<<"number of equivalence classes of proto-synergies "<<xeqc<<endl;
  cout<<"number of equivalence non separables proto-synergies "<<xnsp<<endl;
}

dynamic_bitset<> crn::sp2p(const dynamic_bitset<> sp){
  
  dynamic_bitset<> p(sp_p.size());
  for(int i=0;i<sp_p.size();i++) if(sp_p[i].is_subset_of(sp)) p.set(i,1);
  return(p);
}


dynamic_bitset<> crn::p2sp(const dynamic_bitset<> p){
  
  dynamic_bitset<> sp(lsp);
  for(int i=p.find_first();i!=dynamic_bitset<>::npos;i=p.find_next(i)) sp|=sp_p[i];
  return(sp);
}

void crn::p2b(dynamic_bitset <> &p){
  
  for(int i=0; i<r_p.size();i++){
    if(p[i]!=0) 
      for(int j=p_b[i].find_first();j!=dynamic_bitset<>::npos;j=p_b[i].find_next(j)) 
        if(i!=j) p.set(j,0);
  }
  
}

bool crn::is_sep(const dynamic_bitset<> & sp){
  
  dynamic_bitset<> ib = sp2p(sp);
  p2b(ib);
  
  if(ib.count()<3) return 1;
  
  vector<int> idx = bitget_idx(ib);
  dynamic_bitset<> pb(ib.size());
  int n=ib.count();

  for(int o=1; o<=floor(n/2.0);o++) //o: number of substracted bassics (iterative deepning)
  {

    for(int i=ib.find_first();i!=dynamic_bitset<>::npos;i=ib.find_next(i)){ //i must respect the condition of recusivity end 
      pb.set(i,1);
      if(r_is_sep(o,1,i,ib,pb,sp)) return 1;
      pb.set(i,0);
    }
  }
  return 0;
}

bool crn::r_is_sep(int o, int l, int k, dynamic_bitset<> &ib, dynamic_bitset<> &pb,const dynamic_bitset<> & sp){
//l: nubmer of already substracted  bassics
//k: first substracted canidate (pb index)

  
  if(l==o){
    dynamic_bitset<> csp(sp.size());
    
    csp=p2sp(pb);
    closure(csp);
    if(csp==sp) return 0;
    
    
    pb=pb.flip()&ib;
    csp=p2sp(pb);
    closure(csp);
    
    pb=pb.flip()&ib;
    if(csp==sp) return 0;

    return 1;
  }
  else{
    //tenemos que agregar basicos para verificar separabilidad

    for(int i=ib.find_next(k);i!=dynamic_bitset<>::npos;i=ib.find_next(i)){ 
      
      pb.set(i,1);
      if(r_is_sep(o,l+1,i,ib,pb,sp)) return 1;
      pb.set(i,0);
    }
  }
  return 0;
}

void crn::cont_sep(bool conn=0){
  
  n_cont_clos_p=vector<int>(clos_p.size());
  n_sep_clos_p=vector<int>(clos_p.size());
  ssm=dynamic_bitset<>(clos_p.size());
  rsm=dynamic_bitset<>(clos_p.size());
  
  cout<<"rsm: "<<rsm<<endl;
  cout<<"ssm: "<<ssm<<endl;
  
  if(conn) n_uncon_clos_p= vector<int>(clos_p.size());
  
  for(int i=0;i<clos_p.size();i++){
    
    rsm[i]=is_rsm(clos_sp[i]);
    if(rsm[i]) ssm[i]=is_ssm(clos_sp[i]);

    for(int j=i+1;j<clos_p.size();j++){
      
      if(clos_p[i].is_proper_subset_of(clos_p[j])){
        n_cont_clos_p[j]++;
        
        dynamic_bitset<> b=p2sp(clos_p[j]-clos_p[i]);
        
        if(b!=clos_sp[j]){ 
          n_sep_clos_p[j]++;
          if(conn && !b.intersects(clos_sp[i])) n_uncon_clos_p[j]++;
        }
      }
    }
  }
cout<<"rsm: "<<rsm<<endl;
cout<<"ssm: "<<ssm<<endl;
}
