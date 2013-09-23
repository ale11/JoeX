#ifndef ADT_H
#define ADT_H

//#include <iostream>
#include <vector>
#include <algorithm>

// needed this on uP...
#include <functional>

//#define ADT_LIST_MAX 8192
#define ADT_LIST_MAX 1000000

// ------------------------------------------------------
//note: these where included to satisfy clang mac compiler
using std::max;
using std::min;
// ------------------------------------------------------

template <class T>
class AdtIndexCompare : public std::binary_function<T, T, bool> {
 private:
  const T (*bb)[3];
  int id;                
 public: 
  AdtIndexCompare(const T (*bb)[3],int id) {
    this->bb = bb;
    this->id = id;
  } 

  bool operator()(int index0, int index1) { 
    return( bb[index0][id] < bb[index1][id] );
  }
};

template <class T>
class Adt {

 private:

  typedef struct {
    T bbmin[3];
    T bbmax[3];
    int i_f,i_l;
  } Leafdata;

  int n,nstack,nleaves;
  int * stack;
  Leafdata * leafdata;
  std::vector<int> index;

 public:

  Adt(int n,const T (*bbmin)[3],const T (*bbmax)[3]) {

    if (n < 0) {
      std::cerr << "Error: cannot initialize Adt with n < 0: " << n << std::endl;
      throw(-1);
    }

    initAdt(n,bbmin,bbmax);

  }

  ~Adt() {
    if (stack != NULL) free(stack);
    if (leafdata != NULL) free(leafdata);
  }

  void initAdt(int n,const T (*bbmin)[3],const T (*bbmax)[3]);
  void buildListForBBox(int &n_bbox_list,int * bbox_list,const T bbmin[3],const T bbmax[3]);
  void buildListForPoint(int &n_bbox_list,int * bbox_list,const T point[3]);

};

template <class T>
void Adt<T>::initAdt(int n,const T (*bbmin)[3],const T (*bbmax)[3]) {

  this->n = n;
  stack = NULL;
  leafdata = NULL;
  if (n == 0) return;

  // 1. figure out sizes...
  int i = 1; 
  nstack = 1;
  nleaves = 1;
  while (nleaves < n) {
    // nleaves needs the largest power of 2 greater than n...
    nleaves *= 2;
    // nstack needs 1 + 2 + 3 + 4 + 5 + etc...
    i += 1;
    nstack = nstack + i;
  }
  // and one final multiplier to allow all leaves to
  // contain just one member...
  nleaves *= 2;

  //cout << "Adt n = " << n << " nstack = " << nstack << " nleaves = " << nleaves << std::endl;

  // 2. allocate...
  if ( (stack = (int*)malloc(nstack*sizeof(int))) == NULL ) {
    std::cerr << "Error: memory problem in malloc stack." << std::endl;
    throw(-1);
  }

  if ( (leafdata = (Leafdata*)malloc(nleaves*sizeof(Leafdata))) == NULL ) {
    std::cerr << "Error: memory problem in malloc leafdata." << std::endl;
    throw(-1);
  }

  // add bounding boxes in their passed order...
  index.resize(n);
  for (i = 0; i < n; i++) index[i] = i;

  // set the first leaf's data...
  int ileaf = 0;
  leafdata[ileaf].i_f = 0;
  leafdata[ileaf].i_l = n-1;
        
  int j;
  for (j = 0; j < 3; j++) {
    leafdata[ileaf].bbmin[j] = bbmin[0][j];
    leafdata[ileaf].bbmax[j] = bbmax[0][j];
  }
  for (i = 1; i < n; i++) {
    for (j = 0; j < 3; j++) {
      leafdata[ileaf].bbmin[j] = min(leafdata[ileaf].bbmin[j],bbmin[i][j]);
      leafdata[ileaf].bbmax[j] = max(leafdata[ileaf].bbmax[j],bbmax[i][j]);
    }
  }

  /*
    std::cout << "Adt bbox = " << leafdata[ileaf].bbmin[0] << " < i < " << leafdata[ileaf].bbmax[0] << " , " <<
    leafdata[ileaf].bbmin[1] << " < j < " << leafdata[ileaf].bbmax[1] << " , " <<
    leafdata[ileaf].bbmin[2] << " < k < " << leafdata[ileaf].bbmax[2] << std::endl;
  */

  if (n > 1) {

    // put it on the stack and split recursively...
    stack[0] = ileaf;
    int istack = 1;

    while (istack > 0) {

      // pop the next leaf off the stack...
      ileaf = stack[--istack];

      // get the index range...
      int i_f = leafdata[ileaf].i_f;
      int i_l = leafdata[ileaf].i_l;

      //std::cout << "popping leaf " << ileaf << " i_f/i_l = " << i_f << " " << i_l << std::endl;

      // figure out which direction we should split in -
      // the direction with the largest range...
      int id = 0;
      const T (*bb)[3] = bbmin;
      T delta_max;
      for (j = 0; j < 3; j++) {

        // in the case of bbmin, we store the min in the
        // leaf data, so we only need to compute the max...
        T this_min = leafdata[ileaf].bbmin[j];
        T this_max = bbmin[index[i_f]][j];
        for (i = i_f+1; i <= i_l; i++) this_max = max(this_max,bbmin[index[i]][j]);
        T this_delta = this_max - this_min;
        //std::cout << "bbmin j = " << j << " this_delta = " << this_delta << std::endl;
        if (j == 0) {
          delta_max = this_delta;
        }
        else if (this_delta > delta_max) {
          id = j;
          bb = bbmin;
          delta_max = this_delta;
        }

        // in the case of bbmax, we store the max in the
        // leaf data, so we only need to compute the min...
        this_max = leafdata[ileaf].bbmax[j];
        this_min = bbmax[index[i_f]][j];
        for (i = i_f+1; i <= i_l; i++) this_min = min(this_min,bbmax[index[i]][j]);
        this_delta = this_max - this_min;
        //std::cout << "bbmax j = " << j << " this_delta = " << this_delta << std::endl;
        if (this_delta > delta_max) {
          id = j;
          bb = bbmax;
          delta_max = this_delta;
        }
                        
      }

      //std::cout << "leaf " << ileaf << " splitting in id = " << id << " bb = " << bb << std::endl;

      // sort the index vector from i_f to (and including) i_l according to 
      // bb[index[]][id]...

      /*
        for (i = 0; i < n; i++) std::cout << "before i, index[i] = " << i << " " << 
        index[i] << " " << bb[index[i]][id] << std::endl;
        getchar();
      */

      std::vector<int>::iterator ii_f = index.begin();
      ii_f += i_f;
      std::vector<int>::iterator ii_l = index.begin();
      ii_l += i_l+1; // sort sorts up to but not including second iterator
      std::sort(ii_f,ii_l,AdtIndexCompare<T>(bb,id));
      
      // used to be...
      //std::sort(&(index[i_f]),&(index[i_l+1]),AdtIndexCompare<T>(bb,id));
      
      // then chenged to (to avoid seg faults)...
      //std::sort(&(index[i_f]),&(index[i_l]),AdtIndexCompare<T>(bb,id));
      
      /*
        for (i = 0; i < n; i++) std::cout << "before i, index[i] = " << i << " " << 
        index[i] << " " << bb[index[i]][id] << std::endl;
        getchar();
      */

#ifdef DEBUG
      // check the sort...
      for (i = i_f+1; i <= i_l; i++) {
        if (bb[index[i-1]][id] > bb[index[i]][id]) {
          std::cerr << "Error: sort failed. " << i << " " << i_f << " " << i_l << std::endl;
          throw(-1);
        }
      }
#endif

      // add the children...
      int i_mid = (i_f + i_l)/2;

#ifdef DEBUG
      // check space for 2 children...
      if (2*ileaf+2 >= nleaves) {
        std::cerr << "Error: ileaf out of range." << std::endl;
        throw(-1);
      }
#endif
                        
      // child at 2*ileaf+1 is from i_f to i_mid...
      ileaf = 2*ileaf + 1;
      leafdata[ileaf].i_f = i_f; 
      leafdata[ileaf].i_l = i_mid;
      int this_index = index[i_f];
      for (j = 0; j < 3; j++) {
        leafdata[ileaf].bbmin[j] = bbmin[this_index][j];
        leafdata[ileaf].bbmax[j] = bbmax[this_index][j];
      }
      if (i_mid > i_f) {
        // this leaf still has more than 1 member, so expand
        // bbox and add to stack...
        for (i = i_f+1; i <= i_mid; i++) {
          this_index = index[i];
          for (j = 0; j < 3; j++) {
            leafdata[ileaf].bbmin[j] = min(leafdata[ileaf].bbmin[j],bbmin[this_index][j]);
            leafdata[ileaf].bbmax[j] = max(leafdata[ileaf].bbmax[j],bbmax[this_index][j]);
          }
        }
        if (istack == nstack) {
          std::cerr << "Error: istack out of range." << std::endl;
          throw(-1);
        }
        stack[istack++] = ileaf;
      }

      // child at 2*ileaf+2 is from i_mid+1 to i_l...
      ileaf += 1;
      leafdata[ileaf].i_f = i_mid+1; 
      leafdata[ileaf].i_l = i_l;
      this_index = index[i_mid+1];
      for (j = 0; j < 3; j++) {
        leafdata[ileaf].bbmin[j] = bbmin[this_index][j];
        leafdata[ileaf].bbmax[j] = bbmax[this_index][j];
      }
      if (i_l > i_mid+1) {
        // this leaf still has more than 1 member, so expand
        // bbox and add to stack...
        for (i = i_mid+2; i <= i_l; i++) {
          this_index = index[i];
          for (j = 0; j < 3; j++) {
            leafdata[ileaf].bbmin[j] = min(leafdata[ileaf].bbmin[j],bbmin[this_index][j]);
            leafdata[ileaf].bbmax[j] = max(leafdata[ileaf].bbmax[j],bbmax[this_index][j]);
          }
        }
        if (istack == nstack) {
          std::cerr << "Error: istack out of range." << std::endl;
          throw(-1);
        }
        stack[istack++] = ileaf;
      }

    }

  }

  //for (i = 0; i < n; i++) std::cout << "i, index[i] = " << i << " " << index[i] << std::endl;

  //std::cout << "done Adt()" << std::endl;

}

template <class T>
void Adt<T>::buildListForBBox(int &n_bbox_list,int * bbox_list,const T bbmin[3],const T bbmax[3]) {

  n_bbox_list = 0;

  if (n > 0) {

    // put leaf zero on the stack...
    stack[0] = 0;
    int istack = 1;

    while (istack > 0) {

      // pop the next leaf off the stack...
      int ileaf = stack[--istack];

      // do bbox elimination for the leaf...
      if ( (bbmin[0] > leafdata[ileaf].bbmax[0]) || (bbmax[0] < leafdata[ileaf].bbmin[0]) ||
           (bbmin[1] > leafdata[ileaf].bbmax[1]) || (bbmax[1] < leafdata[ileaf].bbmin[1]) ||
           (bbmin[2] > leafdata[ileaf].bbmax[2]) || (bbmax[2] < leafdata[ileaf].bbmin[2]) )
        continue;

      // this leaf's bbox overlaps. It is the terminal leaf if
      // i_f == i_l...
      if (leafdata[ileaf].i_f == leafdata[ileaf].i_l) {
        // add this index...
        // check...
        if (n_bbox_list >= ADT_LIST_MAX) {
          std::cerr << "Error: increase ADT_LIST_MAX: " << ADT_LIST_MAX << std::endl;
          throw(-1);
        }
        bbox_list[n_bbox_list++] = index[leafdata[ileaf].i_f];
      }
      else {        
        // push both children onto the stack...
        stack[istack++] = 2*ileaf+1;
        stack[istack++] = 2*ileaf+2;
      }

    }

  }

}

template <class T>
void Adt<T>::buildListForPoint(int &n_bbox_list,int * bbox_list,const T point[3]) {

  n_bbox_list = 0;

  if (n > 0) {

    // put leaf zero on the stack...
    stack[0] = 0;
    int istack = 1;

    while (istack > 0) {

      // pop the next leaf off the stack...
      int ileaf = stack[--istack];

      // do bbox elimination for the leaf...
      if ( (point[0] > leafdata[ileaf].bbmax[0]) || (point[0] < leafdata[ileaf].bbmin[0]) ||
           (point[1] > leafdata[ileaf].bbmax[1]) || (point[1] < leafdata[ileaf].bbmin[1]) ||
           (point[2] > leafdata[ileaf].bbmax[2]) || (point[2] < leafdata[ileaf].bbmin[2]) )
        continue;

      // this leaf's bbox overlaps. It is the terminal leaf if
      // i_f == i_l...
      if (leafdata[ileaf].i_f == leafdata[ileaf].i_l) {
        // add this index...
        // check...
        if (n_bbox_list >= ADT_LIST_MAX) {
          std::cerr << "Error: increase ADT_LIST_MAX: " << ADT_LIST_MAX << std::endl;
          throw(-1);
        }
        bbox_list[n_bbox_list++] = index[leafdata[ileaf].i_f];
      }
      else {
        // push both children onto the stack...
        stack[istack++] = 2*ileaf+1;
        stack[istack++] = 2*ileaf+2;
      }

    }

  }

}

#endif
