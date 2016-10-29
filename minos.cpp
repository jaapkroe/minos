#include <vector>
#include <stdio.h>
#include <getopt.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <set>
#include <queue>
#include <unordered_set>
#include <unordered_map>
#include <limits>                 // numeric_limits
#include <ctime>
#include <cmath>                  // round, sqrt
#include <algorithm>
#include "fftw3.h"
#ifdef HASLEMON
#include <lemon/smart_graph.h>
#include <lemon/list_graph.h>
#include <lemon/dijkstra.h>
#include <lemon/bfs.h>
#endif
#include <ctime>                  // clock

#define BONDS 3 // estimate number of bonds to reduce malloc calls (no problem if _v grows larger)
#define MAXLINEW numeric_limits<std::streamsize>::max()
#define DEBUG 0
#define abs(x) (x<0?-x:x)
#define FPREC double
#define IPREC unsigned
#define HASHTYPE unsigned long long 
#define BIG 1e8 // big number for cell size

#ifdef HASLEMON
//#define GR SmartGraph
//#define GR SmartDigraph
#define GR ListDigraph
using namespace lemon;
#endif

using namespace std;

vector<string> elements = {"X","H","He","Li","Be","B","C","N","O","F","Ne","Na","Mg","Al","Si","P","S","Cl","Ar","K","Ca","Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr","Rb","Sr","Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","Sb","Te","I","Xe","Cs","Ba","La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu","Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg","Tl","Pb","Bi","Po","At","Rn","Fr","Ra","Ac","Th","Pa","U","Np","Pu","Am","Cm","Bk","Cf","Es","Fm","Md","No","Lr","Rf","Db","Sg","Bh","Hs","Mt","Ds","Rg"};
vector<FPREC> radii = {1.5, 1.2, 1.4, 3.40, 2.0, 1.7, 1.7, 1.7, 1.52, 1.47, 1.54, 1.36, 1.18, 2.0, 2.1, 1.8, 1.8, 2.27, 1.88, 1.76, 1.37, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 1.63, 1.4, 1.39, 1.07, 2.0, 1.85, 1.9, 1.85, 2.02, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 1.63, 1.72, 1.58, 1.93, 2.17, 2.0, 2.06, 1.98, 2.16, 2.1, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 1.72, 1.66, 1.55, 1.96, 2.02, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 1.86, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0};

struct vertex {
  /* the vertex class holds nodes (atoms) with edges (bonds) */
  //vertex() {size=0;}
  vertex() {neigh.reserve(BONDS);}
  vertex(const vertex& other) { 
    for(unsigned i=0;i<other.neigh.size(); i++) 
      neigh.push_back(other.neigh[i]);
    id = other.id;
  }
  //vector<int> neigh;
  vector<vertex*> neigh;
  //int neigh[BONDS];
  //int size;
  double normal[3];
  int dist, id;
  bool ingraph;
};

struct graph {
  /* the graph class is the main class of this program
  it defines all relevant functions like ring-search */
  graph(const char* fname, double rcut, int verb) : _verb(verb) { // constructor
    _f.open(fname);
    if(!_f) { printf("ERROR opening file %s\n",fname); exit(1); };
    _cell.resize(3);
    _frame=0;
    precompute_r2map(rcut);
  }

  ~graph() {
    _f.close();
  }

  void precompute_r2map(double rcut) {
    /* compute r^2 cutoff between atom types as geometric mean of the radii */
    unsigned n=elements.size(), i, j;
    r2map.reserve(n);
    for(i=0;i<n;i++) {
      r2map[elements[i]].reserve(n);
      for(j=0;j<n;j++) {
        if(rcut>0) r2map[elements[i]][elements[j]] = rcut*rcut;
        else r2map[elements[i]][elements[j]] = radii[i] * radii[j];
      }
    }
    // @modification possible: here one could modify a bond length, e.g.
    //r2map["C"]["C"]=1.67*1.67;        // rcut = 1.67
    //r2map["B"]["B"]=pow(1.45,2);      // rcut = 1.45
    //r2map["B"]["N"]=2.50;             // rcut = 1.58
    //r2map["N"]["B"]=r2map["B"]["N"];  // symmetrize
  }

  int next() {               // read next frame
#ifdef HASLEMON
    _g.clear();
    _nodes.clear();
#endif
    _v.clear();
    _typeset.clear();
    _r2max=0;
    if(loadxyz()) return 1;  // check for EOF (no more frames)
    _frame++;
    connect();
    return 0;
  }

  int loadxyz() {
    string line;
    if(!(_f >> _size)) return 1;                    // read num of atoms (error=EOF!)
    _f.ignore(MAXLINEW, '\n');                      // go to end-of-line
    _types.clear();                                 // clear vectors
    _pos.clear();
    _typesnum.clear();
#ifdef HASLEMON
    //_g.reserveEdge(BONDS*_size);
    _g.reserveNode(_size);
    _g.reserveArc(BONDS*_size);                     // expect 3 arcs per node
#endif
    _v.resize(_size);
    getline(_f,line);                               // second header cell (treat with care)
    stringstream ss(line);
    for(unsigned i=0;i<3;i++) {
      if(!(ss>>_cell[i])) {                         // if no cell was found we set it to very big
        _cell[0]=BIG; _cell[1]=BIG; _cell[2]=BIG;
        break;
      } else { _cell[i] = abs(_cell[i]);
        if(_cell[i]>BIG) {
          fprintf(stderr, "WARNING : cell size (%g) exceeding maximum ",_cell[i]); 
          cerr << BIG << endl; 
          _cell[i] = BIG;
        }
      }
    }
    if(_verb) fprintf(stderr,"CELL = { %20.16g, %20.16g, %20.16g }\n",_cell[0],_cell[1],_cell[2]);
    unsigned ntypes=0, i;
    unordered_map<string, int> map_types_nums;      // map types to a number
    for(i=0;i<_size;i++) {                          // be efficient ...
      getline(_f,line);
      /* fscanf method */
      double x,y,z;
      char t[8];
      sscanf(line.c_str(),"%s%lf%lf%lf",t,&x,&y,&z);
      _types.push_back(t);
      _pos.push_back(x);
      _pos.push_back(y);
      _pos.push_back(z);
      /* fscanf method */

      /* pointer method *
      size_t start = line.find_first_not_of(" \t"); // ignore initial whitespace
      size_t pos = line.find(" ", start);           // ..from there search next space
      char* c = &*(line.begin()+pos+1);             // char*, this pointer will updated to move trough the line
      _pos.push_back(strtod(c+start, &c));          // pointer c is updated in strtod
      _pos.push_back(strtod(c, &c));                // fscanf is ~1.5x slower, ifstream >> num is ~3x slower
      _pos.push_back(strtod(c, NULL));              // we don't care after the last item
      _types.push_back(line.substr(start,pos-start));
      // pointer method */
      if(_verb>1) fprintf(stderr,"x[%5d] = [ % 20.16g, % 20.16g, % 20.16g ]\n",i,_pos[3*i],_pos[3*i+1],_pos[3*i+2]);
      /* determine types here to be more efficient later */
      auto res = _typeset.insert(_types[i]);        // returns a pair of <iter,bool>
      if(get<1>(res)) {                             // successful insert,
        _typesnum.push_back(ntypes);                // set numeric types
        map_types_nums[_types[i]]=ntypes++;         // map: types <-> numeric types
      } else {                                      // no successful insert, key exists
        _typesnum.push_back(map_types_nums[_types[i]]);
      }
#ifdef HASLEMON
      /* add one more node to lemon graph */ 
      GR::Node n = _g.addNode();
      _nodes.push_back(n);
#endif
      _v[i].id=i;
    }

    r2mapnum.clear();
    for(unsigned i=0;i<=ntypes;i++) r2mapnum.push_back(vector<FPREC>(ntypes+1));
    for(auto it=_typeset.begin(); it!=_typeset.end(); ++it ) {
      int i = map_types_nums[*it];
      if(find(elements.begin(),elements.end(),*it)==elements.end()) {
        printf("ERROR : element (%s) unkown\n",it->c_str()); return 2;
      }
      for(auto jt=_typeset.begin(); jt!=_typeset.end(); ++jt ) {
        int j = map_types_nums[*jt];
        FPREC r2 = r2map[*it][*jt];
        if(r2>_r2max)_r2max=r2;
        r2mapnum[i][j] = r2;
        if(_verb) fprintf(stderr, "cutoff radius r[%s,%s]=%.2f\n",(*it).c_str(),jt->c_str(),sqrt(r2));
      }
    }
    return 0;
  }

  
  void king() { 
    //vector<vertex> v(_v); // make a copy of vertex-list, copy-constructor is automatic ... OR IS IT?
    //vector<vertex> v;
    //for(int i=0;i<_v.size();i++) { v.push_back(_v[i]); }
    for(auto vit=_v.begin(); vit!=_v.end(); _v.erase(vit++)) {        // loop over all nodes (remove when done)
      //if(_verb) printf("BFS: vertex \n");
      auto& n = vit->neigh;
      for(auto vjt=_v.begin(); vjt!=_v.end(); vjt++) vjt->ingraph=true;   // reset visitors
      //for(int i=0; i<_v.size(); i++) _v[i].ingraph=true;
      auto jit=n.begin();
      //printf("INGRAPH %d,%d   %d\n" ,_v[0].ingraph,_v[1].ingraph,(*jit)->ingraph);
      //printf("INGRAPH v id    %d\n" ,_v[0].ingraph,_v[1].ingraph,(*jit)->ingraph);
      //break;
      unsigned s = n.size();
      for(unsigned j=0; j<s; j++) {           // loop over all neighbors
        auto save = *jit; 
        n.erase(jit);                           // remove arc i->j
        //printf("v/save in graph ?? = %d, %d\n",vit->ingraph,save->ingraph);
        int success = bfs(*vit,*save,7);        // perform search from v[i]->v[j]
        if(success) printf("cycle found %d<->%d at dist=%d\n",vit->id,save->id,success);
        else printf("no cycle (%d,%d)\n",vit->id,save->id);
        n.push_back(save);                      // reset arc i->j
      }
    }
    return;
  }

  int bfs(vertex& root, vertex& target, int maxdepth) {
    /* breadth-first search : root->target */
    //printf("bfs: %d [nn=%d] -> %d   ingraphs=(%d,%d)\n",root.id,root.neigh.size(),target.id,root.ingraph,target.ingraph);
    //*
    root.dist = 0;
    queue<vertex*> Q; //({root});
    //printf("push root %d dist set to %d\n",root.id,root.dist);
    Q.push(&root);
    printf("\n  root=%d... ",root.id);
    while(!Q.empty()) {
      auto v = Q.front(); 
      printf(" -> %d[%d]",v->id,v->dist);
      Q.pop();
      //printf("popped vertex %d at distance %d [in search for %d]\n",v->id,v->dist,target.id);
      if(v->dist>maxdepth) return 0;
      if(v==&target) return v->dist;
      //if(v.id==target.id) return v.dist;
      //printf("popped vertex %d at distance %d [in search for %d]\n",v->id,v->dist,target.id);
      for(unsigned i=0; i<v->neigh.size(); i++) { //vertex* nit=v.neigh.begin(); nit!=v.neigh.end(); nit++) {
        auto n = v->neigh[i];
        //printf("... neighbor %d , id=%d [ingraph=%d]\n",i,n->id,n->ingraph);
        if(n->ingraph) {
          n->ingraph=false;
          n->dist=v->dist+1;
          Q.push(n);
        }
      }
    }
    //*/
    return 0;
  };

#ifdef HASLEMON
  int lemon_king() { 
    /* king's criterium ring search using the lemon library */
    Bfs<GR> bfs(_g);
    GR::NodeMap<int> dist(_g);
    GR::ArcMap<bool> len(_g,true);
    for(unsigned i=0; i<_size; i++) {
      printf("i=%d\n",i);
      auto s=_nodes[i]; // source
      int j=0;
      for (ListDigraph::ArcIt a(_g); a != INVALID; ++a) {
        printf("j=%d\n",j++);
        _g.reverseArc(a);
        auto t(_g.oppositeNode(s,a));
        //auto t=_nodes[_v[i].neigh[j]->id];
        bfs.run(s,t);
        _g.reverseArc(a);
        //_g.addArc(s,t);
        //s.neigh.push_back(*t);
      }
      _g.erase(s);
    }
  }
#endif


#ifdef HASLEMON
  int lemon_franzblau() {
    if(DEBUG) printf("lemon_franzblau\n");
    Bfs<GR> bfs(_g);
    
    /* ring search using the shortest-path criterion */
    GR::NodeMap<int> dist(_g);
    GR::ArcMap<bool> len(_g);
    bfs.init();
    //for (GR::NodeIt nit(_g); nit!=INVALID; ++nit) ++count;
    //GR::ArcMap<int> length(_g);
    /*
    Dijkstra<GR> dijkstra(_g,length);
    dijkstra.distMap(dist);
    GR::Node s = _nodes[0];
    GR::Node t = _nodes[0];
    //dijkstra.init();
    //dijkstra.addSource(s,10);
    //dijkstra.run();
    //dijkstra.start();
    //dijkstra.run(_nodes[0]);
    for(int i=0;i<_size;i++) {
    //dijkstra.addSource(s,10);
    }
     // dijkstra.run(_nodes[i]);
    //for(int i=0;i<_size;i++)
    //  printf("%i ", dist[_nodes[i]]);
    //printf("\n");
    //Path<GR> p = dijkstra.path(t);
    */
    return 5;
  }
#endif

  int subgraphs() {
    /* detect disconnect subgraphs */
    return 5;
  }
  int shells(int i) {
    /* find neighbor shells from atom index i */
    return 5;
  }

  void connect() {
    /* this function computes the atom neighbors and turns atom positions into graph 
    connections, the algorithm uses verlet lists with maps for sparse connections */
    vector<IPREC> n(3), m(3), m1(3), m2(3), nmax(3);    // cell-indices n,m; limits m1,m2=[-1,0,+1],nmax
    FPREC f[3];                                         // cell size fractions
    typedef unordered_map<long int,vector<int> > mymap; // use hashmap
    HASHTYPE nhash, mhash;                              // hash keywords
    mymap nodemap;                                      // (sparse storage) map between cell id and array of nodes
    mymap::iterator it1, it2;                           // map iterators
    vector<int>::iterator jt1, jt2;                     // vector iterators

    for(unsigned j=0;j<3;j++) {
      /* determine maximum box indices */
      f[j] = 1.0/sqrt(_r2max);
      nmax[j] = (IPREC)(_cell[j]*f[j]);
      if(nmax[j]<3) m1[j]=0; else m1[j]=-1; // for very small cells don't loop too far
      if(nmax[j]<2) m2[j]=0; else m2[j]=1;  // to avoid atoms connected to themselves
      if(_verb) fprintf(stderr,"In direction [%d], #boxes=%d\n",j,nmax[j]);
    }
    for(unsigned i=0;i<_size;i++) {
      /* determine cell for each atom */
      for(unsigned j=0;j<3;j++) {
        double x = fmod(_pos[3*i+j],_cell[j]); // periodic boundary conditions
        if(x<0) x+=_cell[j];
        n[j] = int(x*f[j]);
      }
      hash(n[0],n[1],n[2],nhash);                     // integer hash of boxnumber
      it1 = nodemap.find(nhash);                                // check if key exists
      if(it1==nodemap.end()) nodemap[nhash] = vector<int>(1,i); // add new element to map..
      else nodemap[nhash].push_back(i);                         // ..or push into existing vector
      if(_verb>1) fprintf(stderr,"cell[%3d] = [%5d,%5d,%5d], hash=%lld\n",i,n[0],n[1],n[2],nhash);
    }
    for(it1 = nodemap.begin(); it1 != nodemap.end(); nodemap.erase(it1++)) {
      /* loop over boxes to find neighbors (erase box after each step) */
      nhash = it1->first;
      unhash(n[0],n[1],n[2],nhash); // unhash to get indices
      if(_verb>2) fprintf(stderr,"unhashing %lld  --> (nx,ny,nz) = (%d,%d,%d)\n",nhash,n[0],n[1],n[2]);
      // loop over neighboring boxes
      for(IPREC mx=n[0]+m1[0]; mx<=n[0]+m2[0]; mx++) {
        for(IPREC my=n[1]+m1[1]; my<=n[1]+m2[1]; my++) {
          for(IPREC mz=n[2]+m1[2]; mz<=n[2]+m2[2]; mz++) {
            m[0] = ((mx)%nmax[0]+nmax[0])%nmax[0];                                // positive modulo : mod(x,n) = (x%n+n)%n
            m[1] = ((my)%nmax[1]+nmax[1])%nmax[1];                                // note, otherwise in C(++) -7%3=-1 !!
            m[2] = ((mz)%nmax[2]+nmax[2])%nmax[2];
            //m[0]=mx;
            //m[1]=my;
            //m[2]=mz;
            bool samecell=(n==m);
            if(_verb>2) fprintf(stderr,"cells n<->m = (%u,%u,%u) <-> (%u,%u,%u) [%u]\n",n[0],n[1],n[2],m[0],m[1],m[2],samecell);
            hash(m[0],m[1],m[2],mhash); 
            if(samecell) it2 = it1;
            else it2 = nodemap.find(mhash);
            if(it2 != nodemap.end()) {                                            // if box m exists..
              for(jt1 = it1->second.begin(); jt1 != it1->second.end(); jt1++) {   // jt1-loop (particles in box n)
                for(jt2 = it2->second.begin(); jt2 != it2->second.end(); jt2++) { // jt2-loop (particles in box m)
                  int n1=*jt1, n2=*jt2;
                  if((!samecell&&(n1!=n2)) || (samecell&&(n1<n2))) {              // different cell or n1<n2 (lower half)..
                    FPREC r2 = 0.;                                                // bond distance squared
                    FPREC r2cut = r2mapnum[_typesnum[n1]][_typesnum[n2]];         // r^2 cutoff depending on atom types
                    for(unsigned j=0, j1=3*n1, j2 =3*n2; j<3; j++, j1++, j2++) {
                      FPREC dx = _pos[j1]-_pos[j2];
                      dx -= round(dx/_cell[j])*_cell[j];                          // periodic boundary conditions
                      r2 += dx*dx;
                      if(r2>r2cut) break;
                    }
                    if(r2<r2cut) {                                                // neighbors if smaller than cutoff
#ifdef HASLEMON
                      GR::Arc a = _g.addArc(_nodes[n1],_nodes[n2]);
                      GR::Arc b = _g.addArc(_nodes[n2],_nodes[n1]);
#endif
                      _v[n1].neigh.push_back(&_v[n2]);
                      _v[n2].neigh.push_back(&_v[n1]);
                      //_v[n1].neigh.push_back(n2);
                      //_v[n2].neigh.push_back(n1);
                    }
                    if(_verb>3) fprintf(stderr,"atoms %3d ... %3d : %.2f\n",n1,n2,sqrt(r2));
                  } // n1!=n2
                } // jt2-loop
              } // jt1-loop
            } // m exists?
          } // mx
        }   // my
      }     // mz, neighboring box loops
    } // central box loop
  }

  void hash(IPREC i, IPREC j, IPREC k, HASHTYPE &hash) {
    HASHTYPE ii(i),jj(j),kk(k);
    // packing three 16-bit integers into a single 48-bit (or 64 in reality) integer
    hash = ( ii ) + ( jj << 16 ) + ( kk << 32 ); 
  }

  void unhash(IPREC &i, IPREC &j, IPREC &k, HASHTYPE hash) {
    HASHTYPE len = 65535; // 2^16-1 = 0000 ... 0000 1111 1111 1111 1111 (note preceding zeros in 64-bit binary)
    k = ( hash >> 32 ) & len;  // shift right by 32, then take the last 16 digits in binary
    j = ( hash >> 16 ) & len;  // shift right by 16, then take the last 16 digits in binary
    i = hash & len;            // i is simply the last n binary digits
  }

  void statistics(int lfile) {
    double nnav = 0;
    unsigned i, len;
    vector<int> hist(3);
    for(i=0;i<_size;i++) {
      len = _v[i].neigh.size();
      //int len = countOutArcs(_g,_nodes[i]);
      nnav += len;
      if(len>hist.size()-1) hist.resize(len+1);
      hist[len]++;
    }
    printf("%4d %.4f : ",_frame,nnav/double(_size));
    //printf("%4d %.4f : ",_frame,countArcs(_g)/double(countNodes(_g)));
    for(i=0;i<hist.size();i++) printf("%6d ",hist[i]);
    printf("\n");

    if(lfile) {
      ofstream f("neighbors.dat",ios_base::out);
      f << "NEIGH" << endl;
      for(unsigned i=0;i<_size;i++) {
        f << i << " : ";
        for(unsigned j=0;j<_v[i].neigh.size();j++) {
          f << _v[i].neigh[j]->id << " ";
        }
        f << endl;
      }
      f.close();
    }
  }

  void stat_coordination(int lfile, int nframe) {
    FILE *file=NULL;
    if(lfile) {
      if(nframe==1) file = fopen("coordination.dat","w");
      else file = fopen("coordination.dat","a");
    }
    for(unsigned i=0; i<_size; i++) {
      unsigned n = _v[i].id;
      unsigned nn = _v[i].neigh.size();
      if(!lfile) printf("%-5d %d\n",n,nn);
      else fprintf(file,"%d\n",nn);
    }
    if(lfile) fclose(file);
  }

  void stat_bondlengths(int lfile) {
    FILE *file=NULL;
    if(lfile) file = fopen("bonds.dat","w");
    unsigned i,j,k,n,m;
    for(i=0; i<_size; i++) {
      n = _v[i].id;
      if(lfile<2) printf("%-3d : ",n);
      for(j=0; j<_v[i].neigh.size(); j++) {
        m = _v[i].neigh[j]->id;
        double r2 = 0;
        for(k=0; k<3; k++) {
          double dx = _pos[3*n+k]-_pos[3*m+k];
          dx -= round(dx/_cell[k])*_cell[k]; // periodic boundary conditions
          r2 += dx*dx;
        }
        if(lfile<2) printf (" %3d",m);
        //if(lfile<2) printf (" %3d %7.3f ",m,sqrt(r2));
        if(lfile && n<m) {
          // print average coordinates and bond lengths
          fprintf(file,"%7.3f %7.3f %7.3f    %7.3f\n", 0.5*(_pos[3*n]+_pos[3*m]), 0.5*(_pos[3*n+1]+_pos[3*m+1]), 0.5*(_pos[3*n+2]+_pos[3*m+2]) , sqrt(r2));
        }
      }
      if(lfile<2) printf("\n");
    }
    if(lfile) fclose(file);
  }

  void stat_pyramidalization() {
    unsigned i,j,k;
    for(i=0; i<_size; i++) {
      if(_v[i].neigh.size()!=3) continue;
      int n = _v[i].id;
      printf("%-3d : ",n);
      double p1=0, p2=0, p3=0, p1p3=0, p1p2=0, p2p3=0, theta_p=0;
      double r[3][3]{0};
      for(j=0; j<_v[i].neigh.size(); j++) {
        int m = _v[i].neigh[j]->id;
        for(k=0; k<3; k++) {
          double dx = _pos[3*n+k]-_pos[3*m+k];
          dx -= round(dx/_cell[k])*_cell[k]; // periodic boundary conditions
          r[j][k]=dx;
        }
      }
      for(int k=0; k<3; k++) {
        p1 += r[0][k]*r[0][k];
        p2 += r[1][k]*r[1][k];
        p3 += r[2][k]*r[2][k];
        p1p2 += r[0][k]*r[1][k];
        p1p3 += r[0][k]*r[2][k];
        p2p3 += r[1][k]*r[2][k];
      }
      p1=sqrt(p1); p2=sqrt(p2); p3=sqrt(p3);

      // compute Haddon's pyramidalization angle
      // [R. C. Haddon, J. Am. Chem. Soc., 1997, 119, 1797]
      double dM11 = r[0][0]*(r[1][1]*r[2][2]-r[1][2]*r[2][1]) 
                  - r[1][0]*(r[0][1]*r[2][2]-r[0][2]*r[2][1]) 
                  + r[2][0]*(r[0][1]*r[1][2]-r[0][2]*r[1][1]);
      if(dM11!=0) { // if dM11=0, coordinates are planar
        double dM12 =  p1*(r[1][1]*r[2][2]-r[1][2]*r[2][1]) - p2*(r[0][1]*r[2][2]-r[0][2]*r[2][1]) + p3*(r[0][1]*r[1][2]-r[0][2]*r[1][1]);
        double dM13 =  p1*(r[1][0]*r[2][2]-r[1][2]*r[2][0]) - p2*(r[0][0]*r[2][2]-r[0][2]*r[2][0]) + p3*(r[0][0]*r[1][2]-r[0][2]*r[1][0]);
        double dM14 =  p1*(r[1][0]*r[2][1]-r[1][1]*r[2][0]) - p2*(r[0][0]*r[2][1]-r[0][1]*r[2][0]) + p3*(r[0][0]*r[1][1]-r[0][1]*r[1][0]);

        // coordinates of a sphere going trough r0,r1,r2 and the r_center = (0,0,0)
        double x0 =  0.5*dM12/dM11;
        double y0 = -0.5*dM13/dM11;
        double z0 =  0.5*dM14/dM11;

        double dot_xr0 = x0*r[0][0]+y0*r[0][1]+z0*r[0][2];
        double R=sqrt(x0*x0+y0*y0+z0*z0); // sphere radius
        double theta_sp = acos(dot_xr0/(R*p1))*180./M_PI;
        theta_p = 90 - theta_sp;

        /* this alternative method to compute theta_p does not work
        // [T.P. Radhakrishnan and I. Agranat, Struct. Chem. Vol. 2, 19107 (1990)] 
        // sin(alpha) = p1.p3 / sqrt(p1^2*p3^2)... sin(alpa/2) = sin(0.5*asin(sin(alpha)))
        double a = p1p3/(p1*p3); a=sin(0.5*asin(a));
        double b = p1p2/(p1*p2); b=sin(0.5*asin(b));
        double c = p2p3/(p2*p3); c=sin(0.5*asin(c));
        double sin_theta_sp = 2*a*b*c / sqrt(2*a*a*b*b + 2*a*a*c*c + 2*b*b*c*c - pow(a,4) - pow(b,4) - pow(c,4));
        double theta_p      = asin(sin_theta_sp)*180/M_PI;
        */
      }

      printf ("theta_p = %6.2f\n",theta_p);
    }
  }

  void stat_normal_correlation() {
    /* compute the normal vectors for each atom */
    if(_verb) fprintf(stderr,"computing normal vectors.\n"); 
    compute_normals();
    /* compute forward and backward fourier transforms */
    compute_2d_fft();
  }
    
  void finish_normal_correlation() {
    FILE* f=fopen("correlations.dat","w");
    for(int i=0; i<_nqx; i++) {
      for(int j=0; j<_nqy; j++) {
        int n = i + j*_nqy;
        double correlator = FFT[n][0]*FFT[n][0] + FFT[n][1]*FFT[n][1];
        fprintf (f,"%6d %6d %g\n",i,j,correlator);
      }
      printf("\n");
    }
    fftw_free(FFT);
  }

  void read_grid_indices() {
    char sname[20] = "sample.xyz";
    graph s(sname,-1,_verb);
    if(s.next()) { fprintf(stderr,"ERROR reading %s\n",sname); exit(1); };
    _nqx=round(s.L(0)/s.rmax()); _nqy=round(s.L(1)/s.rmax());

    fprintf(stderr,"read_grid_indices\n");
    char fname[20] = "sample.indices";
    ifstream f;
    string line;
    if(_verb) fprintf(stderr,"nqx, nqy = %d, %d\n",_nqx,_nqy);
    FFT = fftw_alloc_complex(_nqx*_nqy);
    f.open(fname);
    if(!f) { printf("ERROR opening file %s\n",fname); exit(1); };
    while(getline(f,line)) grid_indices.push_back( atoi(line.c_str()) );
  }

  // public access functions
  inline double L(unsigned i) { return (i<_size) ? _cell[i] : 0 ; };
  inline double x(unsigned i, unsigned j) { return (i<_size&&j<3) ? _pos[3*i+j] : 0; };
  inline double xp(unsigned i, unsigned j) { double x=_pos[3*i+j]; return (i<_size&&j<3) ? x-round(x/_cell[j])*_cell[j] : 0; };
  inline double rmax() { return sqrt(_r2max); };
  inline double N() { return _size; };

  void compute_2d_fft() {
    /* 
     * Compute the 2D fourier transforms of the normals
     * and update the correlation function 
     * correlator = <|n(q)|^2> = n(q) * c.c.
     */
    // setup grid sizes
    size_t Nin = 2*_size;
    size_t Nout = _nqx*_nqy;
    unsigned i;
    double *in = fftw_alloc_real(Nin);
    fftw_complex *out = fftw_alloc_complex(Nout);
    fftw_plan p = fftw_plan_dft_r2c_2d(_nqx, _nqy, in, out, FFTW_ESTIMATE);

    // fill input data
    for(i=0; i<Nin; i++) in[i]=0;
    for(i=0; i<_size; i++) {
      // obtain regularized grid point for this atom
      int j = grid_indices[i];
      // store data
      in[2*j]   = _v[i].normal[0];
      in[2*j+1] = _v[i].normal[1];
      //
      in[2*j] = _pos[2*i];
      in[2*j+1] = _pos[2*i+1];
    }

    // run fftw and store results
    fftw_execute(p);
    for(i=0; i<Nout; i++) {
      FFT[i][0] += out[i][0];
      FFT[i][1] += out[i][1];
    }

    fftw_destroy_plan(p); fftw_free(in); fftw_free(out);
  }

  void compute_normals() {
    /* 
     * Normals for atom i are computed as the average 
     * outer-product rij and rik for neigbors j and k of i.
     * An atom needs at least two neighbors to be able to define a normal
     * To ensure the direction is not opposite, every new vector 
     * is projected onto the previously computed normal
     */
    unsigned i,j,k,d;
    double r1[3], r2[3], norm1,norm2;
    double *n;
    for(i=0; i<_size; i++) {
      int numneigh = _v[i].neigh.size();
      double numnorm = 0;
      n = _v[i].normal;
      for(d=0; d<3; d++) n[d]=0;
      if(numneigh<2) continue;
      int ni = _v[i].id;
      if(_verb>1) printf("(NORMAL) %-3d : ",ni);
      for(j=0; j<_v[i].neigh.size(); j++) {
        int nj = _v[i].neigh[j]->id;
        for(d=0,norm1=0; d<3; d++) {
          r1[d] = _pos[3*ni+d]-_pos[3*nj+d];
          r1[d] -= round(r1[d]/_cell[d])*_cell[d]; // periodic boundary conditions
          norm1 += r1[d]*r1[d];
        }
        norm1=sqrt(norm1);
        for(d=0; d<3; d++) r1[d]/=norm1;
        for(k=0; k<j; k++) {
          int nk = _v[i].neigh[k]->id;
          for(d=0,norm2=0; d<3; d++) {
            r2[d] = _pos[3*ni+d]-_pos[3*nk+d];
            r2[d] -= round(r1[d]/_cell[d])*_cell[d]; // periodic boundary conditions
            norm2 += r2[d]*r2[d];
          }
          norm2 = sqrt(norm2);
          for(d=0; d<3; d++) r2[d]/=norm2;
          // compute outer product
          double out_x = r1[1]*r2[2] - r1[2]*r2[1];
          double out_y = r1[2]*r2[0] - r1[0]*r2[2];
          double out_z = r1[0]*r2[1] - r1[1]*r2[0];
          double tmp = out_x*n[0] + out_y*n[1] + out_z*n[2];
          double norm3inv = 1.0/sqrt(out_x*out_x+out_y*out_y+out_z*out_z);
          if(tmp<0) norm3inv = -norm3inv;
          n[0] += (r1[1]*r2[2] - r1[2]*r2[1])*norm3inv;
          n[1] += (r1[2]*r2[0] - r1[0]*r2[2])*norm3inv;
          n[2] += (r1[0]*r2[1] - r1[1]*r2[0])*norm3inv;
          numnorm+=1;
        }
      }
      if(numnorm>0) for(d=0; d<3; d++) n[d] /= numnorm;
      if(_verb>1) printf("%7.3f %7.3f %7.3f\n",n[0],n[1],n[2]);
    }
  }



  inline int frame() { return _frame; }

  private:
    int _nqx, _nqy; 
    fftw_complex *FFT;
    vector<int> grid_indices; // real space mapping of atoms regular grid 
    // for the grid mapping use a reference structure (sample.xyz)
    // e.g.
    // awk 'NR>2{print NR-3,$0}' sample.xyz | sort -k 3,3n -k 4,4n | awk '{print $1}' > sample.indices
    vector<FPREC> _pos, _cell;
    vector<vertex> _v;                // internal storage of vertices
    vector<string> _types;            // atom types
    vector<int> _typesnum;            // numeric version of atom types
    unordered_set<string> _typeset;
    unsigned _size;
    int _frame, _verb;
    ifstream _f;
    FPREC _r2max;                      // maximum r^2[A][B] value for loaded types
    unordered_map<string,unordered_map<string,FPREC> > r2map;  // precomputed r^2[A][B] cutoffs
    vector<vector<FPREC> > r2mapnum;   // precomputed r^2[A][B] cutoffs
#ifdef HASLEMON
    GR _g;
    vector<GR::Node> _nodes;
#endif
};

int main(int argc, char** argv) {
  clock_t start, end;
  char c;

  start = clock();
  string usage="usage: minos <options> file.xyz\n\
  \n\
  options are:\n\
  -a <n>  analysis tool : 1=bond lengths, 2=pyramidalization, 3=coordination, 4=normal correlations\n\
  -f      write analysis data to file\n\
  -c <n>  clusterize defects (defined as rings of size != n)\n\
  -d <n>  ring search depth\n\
  -n <n>  stop after n frames\n\
  -p <n>  print selection (binary code, e.g. 11=1+2+8: 1=neighbors, 2=rings, 4=clusters, 8=chains)\n\
  -r <n>  find rings (criterium: 1=franzblau, 2=king, 3=guttman)\n\
  -R <x>  set manual cutoff distance (default is hardcoded per specied based on vdW radii)\n\
  -s      search connected subgraphs\n\
  -v      increase verbosity level\n\
  -x      print brief neighboring statistics to stdout\n\
  -h      show this help message\n\
  \n\
  Franzblau   : all rings without shortcuts\n\
  King        : all shortest paths between two neighbors\n\
  Guttman     : all shortest paths between root and neighbor\n";

  int  analysis=0;  // analysis tool
  int  clust=0;     // clusterize
  int  depth=7;     // ring search depth
  int  nfrms=0;     // max frame
  int  rings=0;     // ring search (1=franzblau, 2=king, 3=guttman)
  int  print=0;     // print selection 
  bool lstat=false; // print statistics
  int  lfile=0;     // print analysis to file (1=yes, 2=only)
  bool lsubg=false; // subgraph detection
  bool lquiet=false;// be quiet
  int  verbs=0;     // verbosity
  double rcut=-1;   // manual cutoff

  while ((c = getopt(argc, argv, "a:fc:d:n:p:r:R:svqxh")) != -1) {
    switch(c) {
      case 'a': analysis=atoi(optarg); break;
      case 'f': lfile++; break;
      case 'c': clust=atoi(optarg); break;
      case 'd': depth=atoi(optarg); break;
      case 'n': nfrms=atoi(optarg); break;
      case 'p': print=abs(atoi(optarg)); break;
      case 'r': rings=atoi(optarg); break;
      case 'R': rcut=atof(optarg); break;
      case 's': lsubg=!lsubg; break;
      case 'v': verbs++; break;
      case 'q': lquiet=!lquiet; break;
      case 'x': lstat=!lstat; break;
      case 'h': fprintf(stderr,"%s",usage.c_str()); return 0;
      default : fprintf(stderr,"%s",usage.c_str()); return 1;
    }
  }

  if(verbs)
    fprintf(stderr,"# MINOS options: analyze=%d, cluster=%d, depth=%d, print=%d, rings=%d, verbose=%d, stat=%d\n",analysis,clust,depth,print,rings,verbs,lstat);

  if(argc-optind<1) {fprintf(stderr,"ERROR: no input file specified.\n"); return 1;};
  graph g(argv[optind],rcut,verbs);      // initialize graph
  if(analysis==4) g.read_grid_indices();

  while(!g.next()) {   // loop over frames (while not returning an error)
    if(lstat) g.statistics(lfile);   // print statistics
    else cout << "frame " << g.frame() << endl;
#ifdef HASLEMON
    if(rings==1) g.lemon_franzblau();  // find rings
    if(rings==2) g.lemon_king();       // find rings
#else
    if(rings==2) g.king();       // find rings
#endif
    if(lsubg) g.subgraphs();    // find subgraphs
    if(nfrms && g.frame()>=nfrms) break;
    switch(analysis) {
      case 0: break;
      case 1: g.stat_bondlengths(lfile); break;
      case 2: g.stat_pyramidalization(); break;
      case 3: g.stat_coordination(lfile,g.frame()); break;
      case 4: g.stat_normal_correlation(); break;
      default: fprintf(stderr,"ERROR: unimplented analysis method chosen\n"); return 1; break;
    }
  }

  if(analysis==4) g.finish_normal_correlation(); 

  end = clock();
  fprintf(stderr,"CPU time: %.2f sec\n",double(end-start)/CLOCKS_PER_SEC);
}
