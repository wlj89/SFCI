/*

    Sparse Full Configuration Interaction

    L. Wang, Brown University
    
*/

#include "all.h"
#include <cmath>

Hamiltonian H;

unsigned int nThread;
unsigned int blockDim; 
unsigned int threadBcktRatio; 

unsigned int bucketNum;
unsigned int bucketSize;  
unsigned int spaceDim; 
unsigned int bucketSize_bit; 
unsigned int numConnected=20000; 
unsigned int Nd=200000;
unsigned int bucket_unit=60000;


string filename; 
string task;

unsigned int orbitalNum; // number of spin orbitals! 
unsigned int electronNum; 
unsigned int numIteration;
unsigned int numRestart;
DET_TYPE init_hf; 

FLOAT_TYPE target_energy;
FLOAT_TYPE SpMSpVFinalResCut; 
unsigned int HPsiCap;

vector<vector<unsigned int>> nck_list; 

using colvec = dynVec<DET_TYPE,FLOAT_TYPE>; 

void sparse_krylov_solve_eigenval(colvec&,int,int);
void getBasis(colvec&);

void calc_CISD(DET_TYPE, colvec&, int);

FLOAT_TYPE calc_energy(const colvec&);

unsigned long long indexCal(const det_bit&);

void RHF();

void dense_main(); 
void sparse_main(FLOAT_TYPE);

void parse(string key, string val)
{
    if(key == "num_thread")
    {
        nThread = stoi(val);
        cout <<"using input "<< key  << "=" << nThread <<endl ;
    }
    else if(key == "Nd")
    {
        Nd = stoi(val);
        cout <<"using input "<< key  << "=" << Nd <<endl ;
    }
    else if(key == "Nl")
    {
        numIteration = stoi(val);
        cout <<"using input "<< key  << "=" << numIteration <<endl ;
    }   
    else if(key == "num_restart")
    {
        numRestart = stoi(val);  
        cout <<"using input "<< key  << "=" << numRestart <<endl ; 
    }
    else if(key == "num_orb")
    {
        orbitalNum = stoi(val);    
        cout <<"using input "<< key  << "=" << orbitalNum <<endl ; 
        if (orbitalNum > 64)
            throw invalid_argument("orb num larger than 64 is not supported at this point\n");
    }
    else if(key == "num_e")
    {
        electronNum = stoi(val);  
        cout <<"using input "<< key  << "=" << electronNum <<endl ; 
    }
    else if(key == "filename")
    {
        filename = val; 
        cout <<"using input "<< key  << "=" << filename <<endl ; 
    }
    else if(key == "coeff_cut")
    {
        SpMSpVFinalResCut = stod(val); 
        cout <<"using input "<< key  << "=" << SpMSpVFinalResCut <<endl ; 
    }
    else if(key == "bucket_unit")
    {
        bucket_unit = stoi(val); 
        cout <<"using input "<< key  << "=" << bucket_unit <<endl ; 
        //__builtin_popcount() != 2123; 

    }
	else if (key == "init_hf")
	{
		init_hf = stoull(val); 
        cout <<"using input "<< key  << "=" << init_hf <<endl ;
	}
    else if (key == "task")
    {
        task = val;
        cout <<"using input "<< key  << "=" << task <<endl ;
    }
    else if(key == "target_E")
    {
        target_energy = stod(val);
        cout <<"using input "<< key  << "=" << target_energy <<endl ;
    }
}   



void get_paras()
{
    /*
        read parameters from sfci.input
    */
    string raw; 

    ifstream file ("sfci.input");

    if (file.is_open())
    {

        while(getline(file,raw))
        {
            auto pos = raw.find('=');

            string name (raw.begin(),raw.begin() + pos );
            string val (raw.begin() + pos + 1, raw.end()); 
            parse(name,val);
        }

        file.close(); 
    } 
}


int main()
{   
    /*
		You have to correctly set the follwing parameters before running the code. 
		
		Nd: cutoff number 
		blockDim: dimension of spin up/down block 
		nOrb: number of spin orbitals
		nElectron: number of electrons
		bucketSize_bit: bucket size in the SpMSpVM 
		S : the HF determinant in string bit format 	
	*/
	   
    /*
        obtaining parameters
    */
    get_paras();
    
    /*
        loading re-formatted integrals
    */
    load_integrals(filename); 
    cout << "integrals loaded\n";

    unsigned int one = 1; 

    //nThread = 32;       
	
	// range of index. Here only spin-up section is used to calculate index.	
    blockDim = nck(orbitalNum/2,electronNum/2); 
    
	spaceDim = blockDim;	
    threadBcktRatio = 4;    

	bucketSize_bit = ((unsigned int)log2(blockDim))/2; /// 128 

    //avoid expsensive odulus operation bitwise tricks
    bucketSize = one << bucketSize_bit; 

    bucketNum = (spaceDim/bucketSize) + 1; 
	HPsiCap = 60000000;
    numConnected = 10000;   
    SpMSpVFinalResCut = 1e-5; 

	omp_set_num_threads(nThread);
    cout << "\n*****************parameters summary*****************\n";
    cout << "input file:" << filename <<endl;
    cout << "orbital num:" << orbitalNum<<endl;
    cout << "electron num:" << electronNum <<endl;
    cout << "Nd:"<<Nd<<endl;// = 5000000; 
    //cout << "Nd:"<<Nd<<endl;// = 1000000; 
	cout << "core used:" << nThread <<endl;
    cout << "bucket size:" <<bucketSize<<endl;
	cout << "bucket num:" << bucketNum <<endl; 	
    cout << "target energy:" << target_energy<<endl;  
    cout << "SpMSpV final result cut:" << SpMSpVFinalResCut <<endl; 

    // set up the table for nchoosek 
    nck_list.resize(electronNum/2);
    for (int i=0; i< electronNum/2; i++ )
    {
        nck_list[i].resize(orbitalNum/2);
        for (int j =0; j< orbitalNum/2; j++)
        {
            if (j>= i+1)
            {
                nck_list[i][j] = nck(j,i+1); 
            } 
        } 
    }       

	//generateBasis_det(basis);
    //cout <<"Space size:" << basis.size()<<endl;

    spaceDim = blockDim*blockDim; 

    //spaceDim = Nd;
	//spaceDim = 500; 
	//cout <<"FCI space size:" << spaceDim <<endl; 
	cout << "\n****************************************************\n"; 
    
    if (task == "excited")
    {
        /*
            UNDER CONSTRUCTION
        */
        //cout << "Excited state energy calculation starts\n"; 
        //sparse_main(target_energy);
    }
    /*
        
        H2O-6-31g HF: 208911
        c2-cc-pVDZ HF: 202125315
		Ne HF: 12897681423
	    C2-6-31g* HF:202125315 
        C2-6-31g HF:209667 
		N2 cc-pVDZ 26o10e HF: 202125327
		CO : 206208761919
        
		c2 all electron: 4029726735
		F2 all e	   : 857623099392063
	
	*/
	
    /// start lanczos
    else if (task=="ground")
    {
    	colvec st(Nd);
    	    	
    	cout <<"E_hf:" << DiagCal_bit(init_hf) <<endl; 	
    	
    	st.cpyFromPair(init_hf,1.0);	

        auto t1 = system_clock::now(); 
    	//getBasis(st);
    	sparse_krylov_solve_eigenval(st,numIteration,numRestart);

        //sparse_krylov_solve_eigenval_2(st,numIteration,numRestart);
    	
    	auto t2 = system_clock::now();
        std::chrono::duration<double> d1 = t2-t1;
        
        cout << "SFCI time: " << d1.count() <<endl;
    }
    return 0; 	
	
}


/*
!!!!!! Cyber Junkyard !!!!!!
*/

/*
void readPsi(colvec& psi,string filename)
{
    const char *u = &filename[0];
    FILE* fp = fopen(u,"r");

    unsigned int len; 
    fscanf(fp,"%u",&len);
    psi.len = len; 

    DET_TYPE det_temp;
    FLOAT_TYPE amp_temp;
    
    for(unsigned i=0; i<len;i++)
    {
        fscanf(fp,"%llu %lf",&det_temp,&amp_temp);
        psi.basis[i] = det_temp;
        psi.amp[i] = amp_temp; 
    }

    fclose(fp);
}*/
