#include <unistd.h>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <list>
#include <queue>
#include <unordered_set>
#include <unordered_map>
#include <map>
#include <limits>                 // numeric_limits
#include <ctime>                  // clock
#include <cmath>                  // round, sqrt
#include <algorithm>
#include <random>

#define MAXLINEW numeric_limits<std::streamsize>::max()
#define abs(x) (x<0?-x:x)
#define BIG 2e6 // big number for cell size ~ 2^21
using namespace std;

vector<string> elements = {"X","H","He","Li","Be","B","C","N","O","F","Ne","Na","Mg","Al","Si","P","S","Cl","Ar","K","Ca","Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr","Rb","Sr","Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","Sb","Te","I","Xe","Cs","Ba","La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu","Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg","Tl","Pb","Bi","Po","At","Rn","Fr","Ra","Ac","Th","Pa","U","Np","Pu","Am","Cm","Bk","Cf","Es","Fm","Md","No","Lr","Rf","Db","Sg","Bh","Hs","Mt","Ds","Rg"};
vector<double> radii = {1.5, 1.2, 1.4, 3.40, 2.0, 1.7, 1.7, 1.7, 1.52, 1.47, 1.54, 1.36, 1.18, 2.0, 2.1, 1.8, 1.8, 2.27, 1.88, 1.76, 1.37, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 1.63, 1.4, 1.39, 1.07, 2.0, 1.85, 1.9, 1.85, 2.02, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 1.63, 1.72, 1.58, 1.93, 2.17, 2.0, 2.06, 1.98, 2.16, 2.1, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 1.72, 1.66, 1.55, 1.96, 2.02, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 1.86, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0};

struct vertex {
  /* the vertex class holds nodes (atoms) with edges (bonds) */
  //vertex() { neigh.reserve(); } // constructor
  vector<vertex*> neigh;
  string type;
  vector<double> X;
  int id, typenum;
  bool ingraph = true;
};

struct graph {
  /* the graph class defines the ensemble of vertices and functions that can be applied to them */
  graph(const char* fname, double rcut, int verb) : _verb(verb) { // constructor
    _f.open(fname);
    if(!_f) { fprintf(stderr,"ERROR opening file %s\n",fname); exit(1); };
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
    for(i=0; i<n; i++) {
      r2map[elements[i]].reserve(n);
      for(j=0; j<n; j++) {
        if(rcut>0) r2map[elements[i]][elements[j]] = rcut*rcut;
        else r2map[elements[i]][elements[j]] = radii[i] * radii[j];
      }
    }
    // bond length modification possible here, e.g.  
    // r2map["C"]["C"]=1.67*1.67;       // rcut = 1.67
  }

  int next() {               // read next frame
    _v.clear();
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
    _v.clear();
    getline(_f,line);                               // second header cell (treat with care)
    stringstream ss(line);
    for(unsigned i=0; i<3; i++) {
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
    if(_verb>1) fprintf(stderr,"CELL = { %20.16g, %20.16g, %20.16g }\n",_cell[0],_cell[1],_cell[2]);
    unsigned ntypes=0, i;
    unordered_map<string, int> map_types_nums;      // map types to a number
    unordered_set<string> typeset;
    char t[8];
    double x,y,z;
    for(i=0; i<_size; i++) {                        // read coordinates
      getline(_f,line);
      vertex v;
      sscanf(line.c_str(),"%s%lf%lf%lf",t,&x,&y,&z);
      v.X.push_back(x);
      v.X.push_back(y);
      v.X.push_back(z);
      v.type = t;
      v.id = i;
      /* determine unique types */
      auto res = typeset.insert(t);        // returns a pair of <iter,bool>
      if(get<1>(res)) {                     // successful insert,
        map_types_nums[t]=ntypes;           // map: types <-> numeric types
        ntypes++;
      }
      v.typenum = map_types_nums[t];
      _v.push_back(v);
      if(_verb>2) fprintf(stderr,"x[%5d] = [ % 20.16g % 20.16g % 20.16g ]\n",i,x,y,z);
    }

    r2mapnum.clear();
    for(unsigned i=0;i<=ntypes;i++) r2mapnum.push_back(vector<double>(ntypes+1));
    for(auto it=typeset.begin(); it!=typeset.end(); ++it ) {
      int i = map_types_nums[*it];
      if(find(elements.begin(),elements.end(),*it)==elements.end()) {
        printf("ERROR : element (%s) unkown\n",it->c_str()); return 2;
      }
      for(auto jt=typeset.begin(); jt!=typeset.end(); ++jt ) {
        int j = map_types_nums[*jt];
        double r2 = r2map[*it][*jt];
        if(r2>_r2max)_r2max=r2;
        r2mapnum[i][j] = r2;
        if(_verb) fprintf(stderr, "cutoff radius r[%s,%s]=%.2f\n",(*it).c_str(),jt->c_str(),sqrt(r2));
      }
    }
    return 0;
  }

  void connect() {
    /* this function computes the atom neighbors and turns atom positions into graph
    connections, the algorithm uses verlet lists with maps for sparse connections */
    vector<int> n(3), m(3), m1(3), m2(3), nbox(3);                  // cell-indices n,m; limits m1,m2=[-1,0,+1],nbox
    double f[3];                                                     // cell size fractions
    typedef unordered_map<unsigned long long,vector<int> > hashmap; // hashmap
    unsigned long long nhash, mhash;                                // 64-bit hash keywords
    hashmap nodemap;                                                // (sparse storage) map between cell id and array of nodes
    hashmap::iterator it1, it2;                                     // map iterators
    vector<int>::iterator jt1, jt2;                                 // vector iterators

    for(unsigned j=0;j<3;j++) {     // determine maximum box indices
      double eps = 0.01;            // numerical safety margin on box size
      f[j] = 1.0/sqrt(_r2max+eps);
      unsigned long long nbox_long = (unsigned long long)(_cell[j]*f[j]);
      unsigned long long nbox_limit = 2097151; // 2^21-1;
      if(nbox_long > nbox_limit || nbox_long > numeric_limits<unsigned>::max()) {
        fprintf(stderr,"WARNING : cell dimension too large for cutoff for efficient hashing :\n");
        fprintf(stderr,"          nbox=%llu, hash limit=%llu, integer limit=%u\n",nbox_long,nbox_limit,numeric_limits<unsigned>::max());
        unsigned nbox_max = numeric_limits<unsigned>::max() < nbox_limit ? numeric_limits<unsigned>::max()-1 : nbox_limit-1;
        f[j] = nbox_max / _cell[j];
        nbox_long = (unsigned long long)(_cell[j]*f[j]);
      }
      nbox[j] = (int)(nbox_long);
      if(nbox[j]<1) nbox[j]=1;
      if(nbox[j]<3) m1[j]=0; else m1[j]=-1; // for very small cells don't loop too far
      if(nbox[j]<2) m2[j]=0; else m2[j]=1;  // to avoid atoms connected to themselves
      if(_verb) fprintf(stderr,"In direction [%d], #boxes=%d\n",j,nbox[j]);
    }
    if(_verb>1) fprintf(stderr,"loop limits : x[%d,%d] ; y[%d,%d] ; z[%d,%d]\n",m1[0],m2[0],m1[1],m2[1],m1[2],m2[2]);
    for(unsigned i=0;i<_v.size();i++) {            // determine cell for each atom
      for(unsigned j=0;j<3;j++) {
        double x = fmod(_v[i].X[j],_cell[j]); // periodic boundary conditions
        if(x<0) x+=_cell[j];
        n[j] = int(x*f[j]);
        if(n[j]==nbox[j] && nbox[j]>0) n[j]--; // merge the last box with the one before it to avoid a tiny last box
      }
      hash(n[0],n[1],n[2],nhash);              // integer hash of boxnumber
      if(_verb>3) {
        unhash(m[0],m[1],m[2],nhash); // unhash to get indices
        fprintf(stderr,"hashing cell of atom %4d : (%7u,%7u,%7u) HASH-> %20llu UNHASH-> (%7u,%7u,%7u)\n",i,n[0],n[1],n[2],nhash,m[0],m[1],m[2]);
      }
      it1 = nodemap.find(nhash);                                // check if key exists
      if(it1==nodemap.end()) nodemap[nhash] = vector<int>(1,i); // add new element to map..
      else nodemap[nhash].push_back(i);                         // ..or push into existing vector
    }

    for(it1 = nodemap.begin(); it1 != nodemap.end(); nodemap.erase(it1++)) {
      // loop over boxes to find neighbors (erase box after each step)
      nhash = it1->first;
      unhash(n[0],n[1],n[2],nhash); // unhash to get indices
      // loop over neighboring boxes
      for(int mx=n[0]+m1[0]; mx<=n[0]+m2[0]; mx++) {
        for(int my=n[1]+m1[1]; my<=n[1]+m2[1]; my++) {
          for(int mz=n[2]+m1[2]; mz<=n[2]+m2[2]; mz++) {
            m[0] = ((mx)%nbox[0]+nbox[0])%nbox[0];                                // positive modulo : mod(x,n) = (x%n+n)%n
            m[1] = ((my)%nbox[1]+nbox[1])%nbox[1];                                // note, otherwise in C(++) -7%3=-1 !!
            m[2] = ((mz)%nbox[2]+nbox[2])%nbox[2];
            bool samecell=(n==m);
            hash(m[0],m[1],m[2],mhash);
            if(samecell) it2 = it1;
            else it2 = nodemap.find(mhash);
            if(it2 != nodemap.end()) {                                            // if box m exists..
              if(_verb>3) fprintf(stderr,"coupling cells (%7u,%7u,%7u) <-> (%7u,%7u,%7u) [samecell=%u] (#atoms in cell %lu and %lu)\n",n[0],n[1],n[2],m[0],m[1],m[2],samecell,it1->second.size(),it2->second.size());
              for(jt1 = it1->second.begin(); jt1 != it1->second.end(); jt1++) {   // jt1-loop (particles in box n)
                for(jt2 = it2->second.begin(); jt2 != it2->second.end(); jt2++) { // jt2-loop (particles in box m)
                  int n1=*jt1, n2=*jt2;
                  if((!samecell&&(n1!=n2)) || (samecell&&(n1<n2))) {              // different cell or n1<n2 (lower half)..
                    double r2cut = r2mapnum[_v[n1].typenum][_v[n2].typenum];         // r^2 cutoff depending on atom types
                    double r2 = 0.0;
                    for(unsigned d=0; d<3; d++) {
                      double dx = _v[n1].X[d] - _v[n2].X[d];
                      dx -= round(dx/_cell[d])*_cell[d];                          // periodic boundary conditions
                      r2 += dx*dx;
                      if(r2>r2cut) break;
                    }
                    if(r2<r2cut) {                                                // neighbors if smaller than cutoff
                      _v[n1].neigh.push_back(&_v[n2]);
                      _v[n2].neigh.push_back(&_v[n1]);
                    }
                    if(_verb>4) fprintf(stderr,"atoms %3d ... %3d : r=%.2f, rcut=%.2f\n",n1,n2,sqrt(r2),sqrt(r2cut));
                  } // n1!=n2
                } // jt2-loop
              } // jt1-loop
            } // m exists?
          } // mx
        }   // my
      }     // mz, neighboring box loops
    } // central box loop
  }

  void hash(int i, int j, int k, unsigned long long &hash) {
    unsigned long long ii(i),jj(j),kk(k);
    // packing three 21-bit integers into a 63 bits (a 64-bit integer)
    hash = ( ii ) + ( jj << 21 ) + ( kk << 42 );
  }

  void unhash(int &i, int &j, int &k, unsigned long long hash) {
    unsigned long long len = 2097151; // 2^21-1 = 1 1111 1111 1111 1111 1111  (twenty-one 'ones')
    k = ( hash >> 42 ) & len;  // shift right by 42, then take the last 21 digits in binary
    j = ( hash >> 21 ) & len;  // shift right by 21, then take the last 21 digits in binary
    i = hash & len;            // i is simply the last n binary digits
  }

  void statistics() {
    double nnav = 0;
    unsigned i, len, n=0;
    vector<int> hist(3);
    for(i=0; i<_v.size(); i++) {
      if(!_v[i].ingraph) continue;
      len = _v[i].neigh.size();
      nnav += len;
      if(len>hist.size()-1) hist.resize(len+1);
      hist[len]++;
      n++;
    }
    printf("%4d %.4f : ",_frame,nnav/n);
    for(i=0; i<hist.size(); i++) printf("%6d ",hist[i]);
    printf("\n");
  }

  void stat_neighbours() {
    for(auto vi = _v.begin(); vi != _v.end(); vi++) {
      if(!vi->ingraph) continue;
      printf("%d :",vi->id);
      for(auto vj_iter=vi->neigh.begin(); vj_iter!=vi->neigh.end(); vj_iter++) 
        printf(" %d",(*vj_iter)->id);
      printf("\n");
    }
  }

  void stat_coordination() {
    for(unsigned i=0; i<_v.size(); i++) {
      if(!_v[i].ingraph) continue;
      unsigned n = _v[i].id;
      unsigned nn = _v[i].neigh.size();
      printf("%-5d %d\n",n,nn);
    }
  }

  void stat_bondlengths() {
    unsigned i,d,n,m;
    for(i=0; i<_v.size(); i++) {
      n = _v[i].id;
      auto vi = _v[i];
      for(auto vj_iter=vi.neigh.begin(); vj_iter!= vi.neigh.end(); vj_iter++) {
        auto vj = *(*vj_iter); // dereference iterator to pointer
        m = vj.id;
        double r2 = 0;
        for(d=0; d<3; d++) {
          double dx = vi.X[d] - vj.X[d]; 
          dx -= round(dx/_cell[d])*_cell[d]; // periodic boundary conditions
          r2 += dx*dx;
        }
        if(_verb>1) printf ("%.3f %5d %5d %5d %5d\n",sqrt(r2),n,m,vi.typenum,vj.typenum);
        else if(_verb==1) printf ("%.3f %5d %5d\n",sqrt(r2),n,m);
        else printf ("%.3f\n",sqrt(r2));
      }
    }
  }

  void stat_pyramidalization() {
    unsigned i,d;
    for(i=0; i<_v.size(); i++) {
      if(_v[i].neigh.size()!=3) continue;
      int n = _v[i].id;
      printf("%-3d : ",n);
      double p1=0, p2=0, p3=0, p1p3=0, p1p2=0, p2p3=0, theta_p=0;
      double r[3][3]{0};
      int j=0;
      for(auto vj=_v[i].neigh.begin(); vj!= _v[i].neigh.end(); vj++, j++) {
        int m = (*vj)->id;
        for(d=0; d<3; d++) {
          double dx = _v[n].X[d] - _v[m].X[d]; 
          dx -= round(dx/_cell[d])*_cell[d]; // periodic boundary conditions
          r[j][d]=dx;
        }
      }
      for(d=0; d<3; d++) {
        p1 += r[0][d]*r[0][d];
        p2 += r[1][d]*r[1][d];
        p3 += r[2][d]*r[2][d];
        p1p2 += r[0][d]*r[1][d];
        p1p3 += r[0][d]*r[2][d];
        p2p3 += r[1][d]*r[2][d];
      }
      p1=sqrt(p1); p2=sqrt(p2); p3=sqrt(p3);

      // compute Haddon's pyramidalization angle
      // [R. C. Haddon, J. Am. Chem. Soc. 119, 1797 (1997)]
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
      }

      printf ("theta_p = %6.2f\n",theta_p);
    }
  }

  void generate_vacancies(int nvacancies, string t) {
    struct timespec ts; // fast seed
    clock_gettime(CLOCK_MONOTONIC, &ts);
    unsigned int seed = ts.tv_nsec; // current time in nanoseconds
    //seed = 998962097;
    if(_verb) fprintf(stderr,"VACANCIES: SEED = %d\n",seed);
    mt19937 generator(seed); // stdlib Mersenne Twister, 
    for(int i=0;i<nvacancies;i++) {
      uniform_int_distribution<int> distribution(0,_v.size()-1);
      int r = distribution(generator);
      bool selected = false;
      int n=0, nmax=50;
      while(!selected) {
        selected = true;
        // remove an atom where each neighbor is at least 3-coordinated 
        // so as not to leave an 1-coordinated atoms behind
        for(auto n_iter=_v[r].neigh.begin(); n_iter!=_v[r].neigh.end(); n_iter++) {
          if((*n_iter)->neigh.size()!=3) { 
            selected = false;
            break;
          }
        }
        if(!_v[r].ingraph) selected = false;
        if(_v[r].neigh.size()<2) selected = false;
        if(_v[r].type!=t) selected=false;
        if(!selected) r = distribution(generator);
        if(n>nmax) { fprintf(stderr,"ERROR selecting atom : exceeded maximum tries (%d)\n",nmax); exit(1); };
        n++;
      }
      if(_verb) fprintf(stderr,"VACANCIES: selected atom %d (%s)\n",r,_v[r].type.c_str());
      // of those the neighbors of , pick two and move them closer, half the bond length
      uniform_int_distribution<int> pick(0,_v[r].neigh.size()-1);
      int j1 = pick(generator);
      int j2 = pick(generator);
      while(j1==j2) j2 = pick(generator);
      auto n1 = _v[r].neigh[j1];
      auto n2 = _v[r].neigh[j2];
      if(_verb>1) fprintf(stderr,"VACANCIES: bond formation between %d (%d) and %d (%d)\n",j1,n1->id,j2,n2->id);
      double fraction = 0.146; // bond fraction for move
      for(int d=0; d<3; d++) {
        double dx = n1->X[d] - n2->X[d];
        dx -= round(dx/_cell[d])*_cell[d]; // periodic boundary conditions
        dx *= fraction;
        if(dx*dx>1) { fprintf(stderr,"ERROR : movement too large (%f) direction %d : x1, x2 = %f %f\n",dx,d,n1->X[d],n2->X[d]); exit(1); };
        n1->X[d] -= dx;
        n2->X[d] += dx;
      }
      // now we "remove" the vertex : 
      // in reality we just clear it from neighbor lists (remove pointers) and set 'ingraph' to false
      // this avoids dirty pointer business
      // if we simply erase _v[r], it will not only call the destructor on this element but also move all 
      // other elements which is (1) slow and (2) confusion, because now further down in memory are shifted 
      for(auto n_iter = _v[r].neigh.begin(); n_iter!=_v[r].neigh.end(); n_iter++) {
        // remove links to this vertex before removing the vertex itself
        auto nn_iter = (*n_iter)->neigh.begin();
        for(; nn_iter!=(*n_iter)->neigh.end(); nn_iter++) {
          if(*nn_iter == &_v[r]) break; 
        }
        (*n_iter)->neigh.erase(nn_iter);
      }
      _v[r].neigh.clear();
      _v[r].ingraph = false;
      _size--;
    }
  }

  void write_xyz() {
    printf("%d\n%.5f %.5f %.5f\n",_size,_cell[0],_cell[1],_cell[2]);
    for(auto vi = _v.begin(); vi!=_v.end(); vi++) { 
      if(vi->ingraph)
        printf("%-2s % 19.10f % 19.10f % 19.10f\n", vi->type.c_str(), vi->X[0], vi->X[1], vi->X[2]);
    }
  }

  // public access functions
  inline int frame() { return _frame; }

  private:
    vector<double> _cell;
    vector<vertex> _v;                // internal storage of vertices
    int _frame, _verb;
    unsigned _size;
    ifstream _f;
    double _r2max;                      // maximum r^2[A][B] value for loaded types
    unordered_map<string,unordered_map<string,float> > r2map;  // precomputed r^2[A][B] cutoffs
    vector<vector<double> > r2mapnum;   // precomputed r^2[A][B] cutoffs
};

int main(int argc, char** argv) {
  clock_t start, end;
  struct timespec spec;
  time_t s_start, s_end;
  long ns_start, ns_end;
  char c;

  start = clock();
  clock_gettime(CLOCK_MONOTONIC, &spec);
  s_start  = spec.tv_sec; ns_start = spec.tv_nsec;
  string usage="usage: minos <options> file.xyz\n\
  \n\
  options are:\n\
  -a <n>  analysis tool : 1=bond lengths, 2=pyramidalization, 3=coordination, 4=neighbours\n\
  -r <x>  set manual cutoff distance (default is hardcoded per specied based on vdW radii)\n\
  -v      increase verbosity level\n\
  -V <n>  generate n vacancies\n\
  -T <s>  set vacancy type\n\
  -w      write xyz output\n\
  -x      print brief neighboring statistics to stdout\n\
  -h      show this help message\n";

  int  analysis=0;  // analysis tool
  bool lstat=false; // print statistics
  bool lquiet=false;// be quiet
  bool lwrite=false;// write xyz
  int  verbs=0;     // verbosity
  int nvacancies=0; // generate n vacancies
  string tvac = "C";   // vacancy type
  double rcut=-1;   // manual cutoff

  while ((c = getopt(argc, argv, "a:r:V:T:vqwxh")) != -1) {
    switch(c) {
      case 'a': analysis=atoi(optarg); break;
      case 'r': rcut=atof(optarg); break;
      case 'v': verbs++; break;
      case 'V': nvacancies=atoi(optarg); break;
      case 'T': tvac=optarg; break;
      case 'q': lquiet=!lquiet; break;
      case 'w': lwrite=!lwrite; break;
      case 'x': lstat=!lstat; break;
      case 'h': fprintf(stderr,"%s",usage.c_str()); return 0;
      default : fprintf(stderr,"%s",usage.c_str()); return 1;
    }
  }

  if(verbs) fprintf(stderr,"# MINOS options: analyze=%d, verbose=%d, stat=%d\n",analysis,verbs,lstat);

  if(argc-optind<1) {fprintf(stderr,"ERROR: no input file specified.\n"); return 1;};
  graph g(argv[optind],rcut,verbs);      // initialize graph
  while(!g.next()) {                     // loop over frames
    if(nvacancies>0) g.generate_vacancies(nvacancies,tvac);
    if(lwrite) g.write_xyz();
    if(lstat) g.statistics();    // print statistics
    else if(!lquiet) cout << "#frame " << g.frame() << endl;
    switch(analysis) {
      case 0: break;
      case 1: g.stat_bondlengths(); break;
      case 2: g.stat_pyramidalization(); break;
      case 3: g.stat_coordination(); break;
      case 4: g.stat_neighbours(); break;
      default: fprintf(stderr,"ERROR: unimplented analysis method chosen\n"); return 1; break;
    }
  }

  end = clock();
  clock_gettime(CLOCK_MONOTONIC, &spec);
  s_end = spec.tv_sec; ns_end = spec.tv_nsec;
  if(!lquiet) fprintf(stderr,"CPU-time: %.4f sec, Total-time: %.4f sec\n",
        double(end-start)/CLOCKS_PER_SEC, (double) (s_end-s_start) + (double)(ns_end-ns_start)/1e9);
}
