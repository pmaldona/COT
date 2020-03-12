#ifndef _bitcod_
#define _bitcod_

#include <iostream>
#include <boost/dynamic_bitset.hpp>
#include <Rcpp.h>
#include <vector>
#include <cmath>

using namespace std;
using namespace boost;

//futura implementación matrices estequimetricas.
//Coded reaction network Class 
class crn{
public:
  //number of species
  int lsp;
  //number of reactions
  int lrn;
  //bitset of reactats an products for each reactions
  
  vector<dynamic_bitset<> > reac; // reactant species for each reaction
  vector<dynamic_bitset<> > prod; // products species for each reaction
  vector<vector <double> > s_reac; //stochiometry of reactants
  vector<vector <double> > s_prod; //stochiometry of products
  
  
  vector<dynamic_bitset<> > sp_p; // species contained in each partition
  vector<dynamic_bitset<> > sp_b; // species contained in the closure of each partition
  vector<dynamic_bitset<> > sp_sp_b; // stoichiometric positive species contained in the closure of each partition
  vector<dynamic_bitset<> > sn_sp_b; // stoichiometric negative species contained in the closure of each partition
  vector<dynamic_bitset<> > r_p; // reactions supported by each partition (equivalence class)
  vector<dynamic_bitset<> > r_b; // reactions supported by each bassic (equivalence class)
  vector<dynamic_bitset<> > rsp_b; // reactants species contained each bassic (equivalence class)  
  vector<dynamic_bitset<> > p_b; // partitions (equivalence classes) contained in each closure
  vector<int> x_r_p; // index of partitions (equivalence class) that contains the reaction
  
  // vector<dynamic_bitset<> > conn; //which of the non contained basics are connected
  vector<dynamic_bitset<> > conn; //which of the non contained basics are connectable
  vector<vector<dynamic_bitset<> > > mgen; //minimal generator of equivalence classes
  
  vector<dynamic_bitset<> > syn; //(porto) synergies, generating bitsets of partitions 
  vector<dynamic_bitset<> > syn_clos; //clousure of synergies (species) 
  //vector<dynamic_bitset<> > syn_eqc; //equivalence classes of synergies contained in ecah syn_eqc
  vector<int> xsyn_eqc; // index of the equivalence classes to which the synergies belong 
  
  vector<dynamic_bitset<> > reac_part; // reactant partitions (proto-synergies (and monergies)) 
  vector<dynamic_bitset<> > prod_part; // products of proto-synergies (and monergies)
  
  int N_clos; //maximum number of reactive closed sets
  vector<dynamic_bitset<> > clos_sp; // reactive closed sets (species)
  vector<dynamic_bitset<> > clos_p; // reactive closed sets (partitions)
  vector<map <dynamic_bitset<>,vector<vector<int> > > > c_clos_p; //conected reactive closed sets (partitions, by level (number of partitions))
  vector<map <dynamic_bitset<>,vector<vector<int> > > > c_sso_p; //conected reactive semi-selfmaitianed closed sets (partitions, by level (number of partitions))
  int cnt_clos; //counter of reactive closed sets
  vector<int> n_clos_sp; // number of species in each close set
  vector<int> n_clos_p; // numer of partition in each close set
  vector<int> n_cont_clos_p; //number of contained closed sets  
  vector<int> n_sep_clos_p; //number of separable contained closed sets    
  vector<int> n_uncon_clos_p; //number of unconnected contained closed sets
  dynamic_bitset<> rsm; // relational self-mantaining reactive closed sets 
  dynamic_bitset<> ssm; // stochiometry self-mantaining reactive closed sets 
  vector<dynamic_bitset<> > ssmc_p; //semi self-mantaining connected sets (contained basics from c_cs)
  vector<dynamic_bitset<> > ssmc_sp; //semi self-mantaining connected sets (contained species, from c_cs)
  vector<dynamic_bitset<> > ssmo_p; //semi self-mantaining connected sets (contained basics, from c_sso)
  vector<dynamic_bitset<> > ssmo_sp; //semi self-mantaining connected sets (contained species, from c_sso)
  
  
  //Cosntructors
  crn(){}
  crn(int,int);
  crn(Rcpp::List, int, bool = false);
  
  // methods
  
  // Clousure of species set, by reference
  void closure(dynamic_bitset<> &);
  
  //One step closure of species set
  void clos_one(const dynamic_bitset<> &, dynamic_bitset<> &);
  
  //relational self-mantaining of a set
  bool is_rsm(const dynamic_bitset<> &);
  
  //stochiometry self-mantaining of a set
  bool is_ssm(const dynamic_bitset<> &);
  
  // equivalence classes generetor (part of the contructor), generates variables clos, eqc and x_r_p.
  void gen_eqc();
  
  //connected non-contained bassics
  dynamic_bitset<> conn_b(const dynamic_bitset <>&);
  
  //non-contained contributing (to sso) bassics
  dynamic_bitset<> contrib_b(const dynamic_bitset <>&);
  
  //whether sets are connected (species)
  bool are_conn(const dynamic_bitset<>, const dynamic_bitset<>);
    
  // generates mgen (minimal generators), (part of the contructor)
  void gen_mgen();
  
  //closed sets generator
  void cs_gen(int N);
  //recursive call for cs_gen (pb possible bassics sp species)
  bool r_cs_gen(dynamic_bitset<> pb, dynamic_bitset<> sp, int o);
  
  //connected closed set generator
  void c_cs_gen();
  
  //connected semi-selfmantained closed set generator
  void c_sso_gen();
  
  //proto-synergies (and nonergies) generator for an given species set (minimal generator)
  void syn_gen(dynamic_bitset<> sp);
  //recursive proto-synergies (and monergies) generator
  void r_syn_gen(dynamic_bitset<> &b, 
                 int o,
                 vector<dynamic_bitset<> > & sspp,
                 vector<int> &xp);
  
  //all sinergies generator (for all minimal generators)
  void all_syn(bool no_bass);
  
  //species to partitions
  dynamic_bitset<> sp2p(const dynamic_bitset<> sp);
  
  //partitions to species
  dynamic_bitset<> p2sp(const dynamic_bitset<> p);
  
  //partitons to basal partitions
  void p2b(dynamic_bitset <> &p);
  
  // separability of set of species
  bool is_sep(const dynamic_bitset<> & sp);
  
  // recursive separability
  bool r_is_sep(int o, int l, int k,  dynamic_bitset<> & ib, dynamic_bitset<> & pb, const dynamic_bitset<> & sp);
  
  // contention reactive contained closed  set, their separability and connectivity
  void cont_sep(bool);
  
};

#endif
