#ifndef ALL_H_INCLUDED
#define ALL_H_INCLUDED

#include <iostream>
#include <fstream>
#include <unordered_map>
#include <unordered_set>
#include <iomanip>
#include <vector>
#include <array>
#include <string>
#include <tuple>
#include <functional>
#include <chrono>
#include <algorithm>
#include <cmath>
#include <ctime>
#include <cstdlib>
#include <Eigen/Eigen>
#include <Eigen/Dense>
#include <Eigen/SparseLU>
#include <unsupported/Eigen/IterativeSolvers>
#include <random>
#include <omp.h>
#include <bitset>
#include <stdexcept>



/*
    Molecule      #e-      #spin-orbs   HF

    H2            2        8
    LiH           2        20
    Li2           2        32
    O2            12       32
    H2O           8        24
    C2H4          12       48
    Cl2           14       32

    UO2           5        28

    UF3BF3        11       28      {1 2 9 10 14 16 18 20 22 24 26}
    UF3BF3-1      8        28      {2 14 16 18 20 22 24 26}

    UF3BF3-2      7        28      {0 1 16 18 22 24 26}
    UF3BF3-2-1    8        28      {0 1 2 3 5 7 12 13}

	Ne aug-cc-pVDZ 10     46
*/


using namespace Eigen;
using namespace std;
using namespace chrono; 


extern unsigned int Nd; 
extern string filename; 
extern unsigned int orbitalNum;
extern unsigned int electronNum; 
extern unsigned int numIteration;
extern unsigned int numRestart;
extern unsigned int bucket_unit;


const unsigned numIterations = 5;   
const unsigned numRestarts = 1000;
const int blockSize = 10000;
const int blockNum = 1;
const int enumNum = 4;
const double bar=0.00001;

//const unsigned int Nd = 400000;

typedef double FLOAT_TYPE; 

//typedef array<unsigned char,electronNum> det;
typedef array<unsigned char,4> twoBody;
typedef array<unsigned char,2> oneBody;
typedef uniform_int_distribution<int> UID;
typedef MatrixXd M;
typedef SparseMatrix<FLOAT_TYPE,RowMajor> SpMat;
typedef Matrix<FLOAT_TYPE, Dynamic, 1> ColVec;
typedef Triplet<FLOAT_TYPE> T;
typedef SelfAdjointEigenSolver<MatrixXd> SAES;
typedef GeneralizedSelfAdjointEigenSolver<MatrixXd> GSAES; 

//typedef unsigned __int128 det_bit; 
typedef unsigned long long det_bit;  
typedef det_bit DET_TYPE; 

/*
class detHash 
{
    public: 
	size_t operator()(const det& rhs) const
	{
		string shadow(rhs.begin(),rhs.end());
		return hash<string>()(shadow);
	}
};
*/

struct Hamiltonian
{
    /// separate dig and off diag operators
    /// use n + orbNum*q (n<q) as the key

    ///Off diagonal part
    vector<vector<twoBody>> OO;
    vector<vector<FLOAT_TYPE>> OO_coeff;

    vector<vector<oneBody>> O;
    vector<vector<FLOAT_TYPE>> O_coeff;

    ///Diagonal part
    vector<twoBody> OO_diag;
    vector<FLOAT_TYPE> OO_coeff_diag;

    vector<oneBody> O_diag;
    vector<FLOAT_TYPE> O_coeff_diag;

    FLOAT_TYPE constant;
};

//typedef unordered_set<det,detHash> UOS;
//typedef unordered_map<det,FLOAT_TYPE,detHash> UOM;

extern Hamiltonian H;
//extern vector<det> basis; 
extern vector<det_bit> basis_bit; 
extern double lambda; 
//extern unordered_map<det,int, detHash> basis_pos; 
extern unsigned int nThread; 
extern unsigned int blockDim; 

extern unsigned int bucketNum;
extern unsigned int bucketSize;
extern unsigned int bucketSize_bit;
extern unsigned int threadBcktRatio; 
extern unsigned int numConnected; 
extern unsigned int HPsiCap;

extern vector<vector<unsigned int>> nck_list; 
extern FLOAT_TYPE SpMSpVFinalResCut;

extern DET_TYPE init_hf;

template <typename DET, typename COEFF>
class dynVec
{   
    /*
        data structure for the wave function 
    */

public:
    unsigned int len;    // actual # of determinant 
    unsigned int cap;    // capacity 

    DET *basis = NULL;
    COEFF *amp = NULL; 

    dynVec()
    {
        len = 0;
        cap = 0; 
    } 

    dynVec(const unsigned int resNum)
    {
        /// constructor 
        len = 0;
        cap = resNum;

        basis = new DET[resNum];
        amp = new COEFF[resNum];

    }      

    void reserveMem(const unsigned int resNum)
    {
        len = 0;
        cap = resNum; 

        basis = new DET[resNum];
        amp = new COEFF[resNum];    
    }

    ~dynVec()
    {
        delete [] basis; 
        delete [] amp; 
    }   
    
    void checkOrder()
    {
        for (int i=1; i<len; i++)
        {
            if (basis[i] < basis[i-1])
            {
                cout << "WARNING: not in ascending order!!!"; 
                break; 
            }

        }
    }

    void cpyFromVec(const dynVec<DET,COEFF>& src)
    {
        #pragma omp simd 
        for(int i=0; i< src.len; ++i)
        {
            basis[len+i] = src.basis[i];
            amp[len+i] = src.amp[i];    
        }
        len += src.len; 
    }
    void cpyFromPair(const DET src_det, const COEFF src_coeff) 
    {   
        basis[len] = src_det;
        amp[len] = src_coeff;
        len ++;   
    }   
    
    void swapIdx(int a, int b)
    {
        swap(basis[a],basis[b]);
        swap(amp[a],amp[b]);
    }

    void pasteFromIdx(int a, int b)
    {
        basis[a] = basis[b];
        amp[a] = amp[b]; 
    }

    void assignValue(int a, DET_TYPE value_basis, FLOAT_TYPE value_amp)
    {
        basis[a] = value_basis;
        amp[a] = value_amp; 
    }

    void setZeroVec()
    {
        len = 0; 
    }

    FLOAT_TYPE norm()
    {
        FLOAT_TYPE n = 0; 

        #pragma omp simd reduction(+ : n)
        for(int i=0; i<len; i++)
            n += amp[i]*amp[i]; 

        return sqrt(n);
    }

    void scalarMtply(FLOAT_TYPE alpha)
    {
        #pragma omp simd
        for(int i=0; i< len; i++)
            amp[i] *= alpha; 
    } 

    void normalize() 
    {
        FLOAT_TYPE n = norm();
        scalarMtply(1.0/n);  
    }
    
    void shrinkInPlace(const unsigned int len)
    {
        void introsort_amp(dynVec<DET_TYPE,FLOAT_TYPE>&,int);
        void introsort(dynVec<DET_TYPE,FLOAT_TYPE>&,int);

        introsort_amp(*this,len);   
    }

    void shrinkFrom(dynVec<DET,COEFF>& src)
    {
        void introsort_amp(dynVec<DET_TYPE,FLOAT_TYPE>&,int);
        void introsort(dynVec<DET_TYPE,FLOAT_TYPE>&,int);

        introsort_amp(src,src.len);

        int st_idx = src.len < cap ? 0 : src.len - cap;
        
        for(unsigned int i = st_idx; i< src.len; i++)
        {
            basis[i - st_idx] = src.basis[i];
            amp[i - st_idx] = src.amp[i];
        } 

        len = src.len - st_idx; 

        introsort(*this,len); 
    }

    void merge()
    {
        /*
            remove determinant duplicacy 
        */
        unsigned int merge_bucket(dynVec<DET_TYPE,FLOAT_TYPE>&, unsigned int); 
        void introsort(dynVec<DET_TYPE,FLOAT_TYPE>&,int);

        introsort(*this,len);

        unsigned int merged_size = merge_bucket(*this,len);

        len = merged_size; 

    }

    FLOAT_TYPE dot(const dynVec<DET,COEFF>& src)
    {   
        
        // num of element is written in len
        if(len ==0 || src.len ==0 ) return 0.0; 

        unsigned int idx_this = 0;
        unsigned int idx_src = 0;   

        FLOAT_TYPE sum = 0.0;   
		//cout << len <<' '<<src.len<<endl;

        while (idx_this < len && idx_src < src.len)
        {
            if (basis[idx_this] < src.basis[idx_src])
            {
                idx_this ++; 
            }
			else 
			{
				if (basis[idx_this] > src.basis[idx_src])
				{
					idx_src ++; 
				}
				else 
				{
					sum += amp[idx_this] * src.amp[idx_src];
					idx_this ++; 
					idx_src ++; 
				}
			}	
        }  
	    //cout << idx_this <<" " <<idx_src<<endl; 	
        return sum; 
    }  

    void addTwoCorrected(const dynVec<DET,COEFF>& w, FLOAT_TYPE a)
    {
        /*
            determinants in w.basis are excluded
        */
        /*
        void introsort_amp(dynVec<DET_TYPE,FLOAT_TYPE>&,int);
        void introsort(dynVec<DET_TYPE,FLOAT_TYPE>&,int);
        unsigned int mergeBuffer(dynVec<DET_TYPE,FLOAT_TYPE>&, unsigned int); 
        */

        dynVec<DET,COEFF> temp; 
        temp.reserveMem(w.len + len); 

        unsigned int idx_this = 0;
        unsigned int idx_w = 0; 

        while(idx_this < len && idx_w < w.len)
        {
            if(basis[idx_this] == w.basis[idx_w])
            {
                idx_this ++;
                idx_w ++; 
            }
            else if (basis[idx_this] > w.basis[idx_w])
            {
                temp.cpyFromPair(basis[idx_this],amp[idx_this]);
                idx_this ++; 
                idx_w ++; 
            }
            else
            {
                temp.cpyFromPair(basis[idx_this],amp[idx_this]);
                idx_this ++; 
            }

        }

        while(idx_this < len)
        {
            temp.cpyFromPair(basis[idx_this],amp[idx_this]);
            idx_this ++; 
        }

        // put things back 

        for(int i=0; i< temp.len; i++)
        {
            basis[i] = temp.basis[i];
            amp[i] = temp.amp[i];
        }
        len = temp.len; 
    }
    
    void addTwo(const dynVec<DET,COEFF>& w, FLOAT_TYPE a)
    {
        /*
            v = v+ a*w , and truncate. 
            No need to sort at the first place 
        */
		//cout << "beep"<<endl;
        void introsort_amp(dynVec<DET_TYPE,FLOAT_TYPE>&,int);
        void introsort(dynVec<DET_TYPE,FLOAT_TYPE>&,int);
        unsigned int merge_bucket(dynVec<DET_TYPE,FLOAT_TYPE>&, unsigned int); 

        dynVec<DET,COEFF> temp; 
        temp.reserveMem(w.len + len);   
		 
        unsigned int idx_this = 0;
        unsigned int idx_w = 0;     
        
        while( idx_this < len && idx_w < w.len)
        {
            /*
                copying two determinant-major sorted arrays
                to form the 3rd that is still ordered
            */  

            if(basis[idx_this] < w.basis[idx_w])
            {
                temp.cpyFromPair(basis[idx_this],amp[idx_this]);
                idx_this ++; 
            }
            else 
            {
                temp.cpyFromPair(w.basis[idx_w],a*w.amp[idx_w]);
                idx_w ++; 
            }
        }
		//cout <<"beep" <<endl;
        while(idx_this < len)
        {   
            temp.cpyFromPair(basis[idx_this],amp[idx_this]); 
            idx_this ++; 
        }

        while(idx_w < w.len)
        {   
            temp.cpyFromPair(w.basis[idx_w],a*w.amp[idx_w]);
            idx_w ++; 
        }
        
		//cout << "asd"<<endl;
        unsigned int merged_len = merge_bucket(temp,temp.len);  
        introsort_amp(temp,merged_len);   

        //taking the top ones in
        unsigned int st = merged_len < cap ? 0 : merged_len - cap; 

        int copy_idx = 0; 
        for (int i = st; i < merged_len; i++)
        {   
            /*
                determinants with 0 amplitude can introdcue serious problem
            */    
            //if (abs(temp.amp[i]) > 1e-14)
            //{
            basis[i-st] = temp.basis[i]; 
            amp[i-st] = temp.amp[i]; 
                //copy_idx ++; 
            //}       
        } 
        len = merged_len - st;  
        //len = copy_idx; 

        introsort(*this,len);  
		//cout << "addTwo basisi sort done\n\n";
    }

    void addThree(const dynVec<DET,COEFF>& x,
                  const dynVec<DET,COEFF>& y,
                  const dynVec<DET,COEFF>& z,
                  FLOAT_TYPE a,
                  FLOAT_TYPE b,
                  FLOAT_TYPE c)
    {
        /*
            v = a*x + b*y + c*z
        */

        void introsort_amp(dynVec<DET_TYPE,FLOAT_TYPE>&,int);
        void introsort(dynVec<DET_TYPE,FLOAT_TYPE>&,int);
        unsigned int merge_bucket(dynVec<DET_TYPE,FLOAT_TYPE>&, unsigned int); 

        dynVec<DET,COEFF> temp;  
        temp.reserveMem(x.len + y.len + z.len); 
        
        unsigned int idx_x = 0;
        unsigned int idx_y = 0;
        unsigned int idx_z = 0;
        
        while (idx_x < x.len && idx_y < y.len && idx_z < z.len)
        {
            DET_TYPE m = min(x.basis[idx_x],min(y.basis[idx_y],z.basis[idx_z]));

            if (m == x.basis[idx_x])
            {
                temp.cpyFromPair(x.basis[idx_x],a*x.amp[idx_x]);
                idx_x++; 
            }
            else if (m == y.basis[idx_y])
            {
                temp.cpyFromPair(y.basis[idx_y],b*y.amp[idx_y]);
                idx_y++;
            }
            else 
            {
                temp.cpyFromPair(z.basis[idx_z],c*z.amp[idx_z]);
                idx_z++;
            }   
        }  
        
        while(idx_x < x.len && idx_y < y.len)
        {
            if(x.basis[idx_x] < y.basis[idx_y])
            {
                temp.cpyFromPair(x.basis[idx_x], a*x.amp[idx_x]);
                idx_x++; 
            }
            else
            {
                temp.cpyFromPair(y.basis[idx_y], b*y.amp[idx_y]);
                idx_y++; 
            }
        }   

        while(idx_x < x.len && idx_z < z.len)
        {
            if(x.basis[idx_x] < z.basis[idx_z])
            {
                temp.cpyFromPair(x.basis[idx_x], a*x.amp[idx_x]);
                idx_x++; 
            }
            else
            {
                temp.cpyFromPair(z.basis[idx_z], c*z.amp[idx_z]);
                idx_z++; 
            }
        }

        while(idx_z < z.len && idx_y < y.len)
        {
            if(y.basis[idx_y] < z.basis[idx_z])
            {   
                temp.cpyFromPair(y.basis[idx_y], b*y.amp[idx_y]);
                idx_y++; 
            }
            else 
            {
                temp.cpyFromPair(z.basis[idx_z], c*z.amp[idx_z]);
                idx_z++; 
            }
        }

        while(idx_x < x.len)
        {
            temp.cpyFromPair(x.basis[idx_x], a*x.amp[idx_x]);
            idx_x++;
        }  

        while(idx_y < y.len)
        {
            temp.cpyFromPair(y.basis[idx_y], b*y.amp[idx_y]);
            idx_y++;
        }  

        while(idx_z < z.len)
        {
            temp.cpyFromPair(z.basis[idx_z], c*z.amp[idx_z]);
            idx_z++;
        }   

        unsigned int merged_len = merge_bucket(temp,temp.len);  
        introsort_amp(temp,merged_len);   
		//cout << "addThree amp sort done\n";
        // updating len of self object 
        unsigned int st = merged_len < Nd ? 0 : merged_len - Nd; 

        int copy_idx = 0; 
        for (int i = st; i < merged_len; i++)
        {   
            /*
                determinants with 0 amplitude can introdcue serious problem
            */    
			if (abs(temp.amp[i]) > 1e-14)
            {
                basis[copy_idx] = temp.basis[i]; 
                amp[copy_idx] = temp.amp[i]; 
                copy_idx ++; 
            }       
        }  
        len = copy_idx;
        
        introsort(*this,len);  
		//cout << "addThree basis sort done\n";
    }

    void naiveCopy(const dynVec<DET,COEFF>& src)
    {
        /*
            Copy without boundary check 
        */
        len = src.len; 

        //#pragma omp simd
        for(unsigned int i=0; i< len; i++)
        {
            basis[i] = src.basis[i];   
            amp[i] = src.amp[i]; 
        }   
    }   

}; 


void load_integrals(string);
//void printDet(const det&);
//void generateBasis_det(vector<det>&);
//void OffDiagGen_det(const det&, vector<det>&,vector<double>&);
//void solve(SpMat&, ColVec&, int);
//void loadFromT(SpMat&, vector<det>&, int);
//void minres(ColVec&, ColVec&, ColVec&, int&, double&); 
void minres_sp(dynVec<DET_TYPE, FLOAT_TYPE>&,
               dynVec<DET_TYPE, FLOAT_TYPE>&,
               dynVec<DET_TYPE, FLOAT_TYPE>&,
               int,FLOAT_TYPE,FLOAT_TYPE); 

void SpMVMtply(ColVec&, ColVec&, double, double); 
//void SpMVMtply_bucket(vector<DET_TYPE>&, vector<DET_TYPE>&, ColVec&, ColVec& , double, double); 
void SpMVMtply_bucket(dynVec<DET_TYPE,FLOAT_TYPE>&, dynVec<DET_TYPE,FLOAT_TYPE>&, double, double,FLOAT_TYPE&); 

//double DiagCal_det(const det&);
double residue(SpMat&, ColVec&); 
double expectation(SpMat&, ColVec&);

FLOAT_TYPE calEnergy(dynVec<DET_TYPE,FLOAT_TYPE>&); 

//det bitToDet(const det_bit);
//det_bit detToBit(const det);
double signCal_bit(const det_bit&, int, int);
//void OffDiagGen_bit(const det_bit, vector<det_bit>&, vector<double>&);
void OffDiagGen_bit(const DET_TYPE, const FLOAT_TYPE, vector<dynVec<DET_TYPE,FLOAT_TYPE>>&, unsigned int);
double DiagCal_bit(const det_bit);

unsigned int signatureCal(const det_bit&);
unsigned int nck(int, int); 
unsigned long long indexCal(const det_bit&);
//int spinNumCal(const det&);




#endif // CLASS_H_INCLUDED
