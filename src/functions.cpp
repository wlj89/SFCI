#include "all.h"

void introsort(dynVec<DET_TYPE,FLOAT_TYPE>&, int); 
void introsort_amp(dynVec<DET_TYPE,FLOAT_TYPE>&, int);

int spinCal_bit(const det_bit);

/*
int spinNumCal(const det& q)
{   
    //
    //    M_s = 0 block only
    //
    int u =0;
    
    for (int i=0; i< q.size(); i++)
    {
        if (q[i] %2 ==0 )
            u ++;
        else 
            u--;
    }
    return u; 
}
*/

/*
void comb_det(int n, int r, int *arr, int sz, vector<det>& basis)
{
    
    //    n = total number
    //    r = number to choose
    //    sz = r
		
    //    define arr before hand
    
    for (int i = n; i >= r; i --)
    {
        arr[r - 1] = i;
        if (r > 1)
            comb_det(i - 1, r - 1, arr, sz, basis);
        else
        {
            det chi {};

            for (int i = 0; i < sz; i ++)
                chi[i] =arr[i]-1;
            
            if (spinNumCal(chi) == 0)
                basis.push_back(chi);
        }
    }
}
*/

void comb_RHF(int n, int r, int *arr, int sz, vector<DET_TYPE>& basis)
{
    /**
        n = total number
        r = number to choose
        sz = r
        
        define arr before hand
    */
    for (int i = n; i >= r; i --)
    {
        arr[r - 1] = i;
        if (r > 1)
            comb_RHF(i - 1, r - 1, arr, sz, basis);
        else
        {
            DET_TYPE one = 1; 
            DET_TYPE result = 0;

            for (int i = 0; i < sz; i ++)
            {
                result += one << (2*(arr[i]-1));
                result += one << (2*(arr[i]-1)+1);
            }
            
            basis.push_back(result);
        }
    }
}

void RHF()
{

    /*
        a stupid restricted Hartree fock 
        Sz = 0 
    */
    vector<DET_TYPE> basis; 

    int n = orbitalNum/2;
    int r = electronNum/2;
    int sz = r;
    int arr[electronNum];   

    comb_RHF(n, r, arr, sz, basis);

    double optimal_diag = 1000.0; 
    DET_TYPE optimal_det; 

    for (int i=0; i< basis.size(); i++ )
    {
        cout << basis[i] <<endl; 

        double diag = DiagCal_bit(basis[i]);

        cout << diag <<endl<<endl; 

        if (diag < optimal_diag)
        {
            optimal_diag = diag;
            optimal_det = basis[i];
        }
    }

    cout << optimal_diag <<endl;
    cout << optimal_det <<endl; 
    cout << DiagCal_bit(optimal_det) << endl ;

}

/*
void generateBasis_det(vector<det>& basis)
{

    int n = orbitalNum;
    int r = electronNum;
    int sz = r;
    int arr[electronNum];

    comb_det(n,r, arr, sz, basis);

}*/

/*
double signCal(det& reff_ord)
{
    int u =0;

    //#pragma omp simd reduction(+:u) collapse(2)
    for(int i=0;i<electronNum-1;i++)
        for (int j=i+1; j<electronNum; j++)
			u = (reff_ord[i] > reff_ord[j]?u+1:u);

    return u%2 ==0 ? 1.0 : -1.0;
}
*/


/*
void OffDiagGen_det(const det& reff, vector<det>& basis_det, vector<double>& amp)
{
    void printDet(const det&);
    
    array<short,orbitalNum> occup;
    for (int i=0;i<orbitalNum; i++)
        occup[i] = -1;
    
    for(int i=0; i< electronNum; i++)
        occup[reff[i]] = i;

    /// two body operators
    for (int i=0; i<electronNum; i++)
        for(int j=i+1;j<electronNum; j++)
        {
            unsigned idx = reff[i]  + reff[j] * orbitalNum;
            unsigned len = H.OO[idx].size();
            
            for(int k = 0; k< len; k++)
            {
                
                det start = reff;

                unsigned char m = H.OO[idx][k][0];
                unsigned char n = H.OO[idx][k][1];
                unsigned char p = H.OO[idx][k][2];
                unsigned char q = H.OO[idx][k][3];

                ///double excitation
                if ( (m!=n) && (n!=p) && (p!=q) && (m!= q) &&
                     (occup[m] ==-1) && (occup[p] ==-1) ) //&& (occup[n] < occup[q]))
                {
                    
                    start[occup[q]] = p;
                    start[occup[n]] = m;

                    //auto start_2 = start;
                    double sign = signCal(start);
                    sort(start.begin(), start.end());

                    basis_det.push_back(start);
                    amp.push_back(sign*H.OO_coeff[idx][k]); 
                    
                }
                
                ///single excitation
                else if( m!=n && n!=p && p==q && q!=m && occup[m] ==-1 )
                {
                    start[occup[n]] = m;

                    //auto start_2 = start;
                    double sign = signCal(start);
                    sort(start.begin(), start.end());

                    basis_det.push_back(start);
                    amp.push_back(sign*H.OO_coeff[idx][k]); 
                }

                else if ( m!=n && n==p && p!=q && q!=m && occup[m] ==-1 )
                {
                    start[occup[q]] = m;

                    //auto start_2 = start;
                    double sign = signCal(start);
                    sort(start.begin(), start.end());

                    basis_det.push_back(start);
                    amp.push_back(-1.0*sign*H.OO_coeff[idx][k]); 

                }

            }

        }
    /// one body operators

    for (int i=0; i<electronNum; i++)
    {
        unsigned idx = reff[i];
        unsigned len = H.O[idx].size();

        for (int j = 0; j< len; j++)
        {
            det start = reff;
            unsigned char m = H.O[idx][j][0];
            unsigned char n = H.O[idx][j][1];
            
            if (occup[m] == -1)
            {
                start[occup[n]] = m;

                double sign = signCal(start);
                sort(start.begin(), start.end());

                basis_det.push_back(start);
                amp.push_back(sign*H.O_coeff[idx][j]); 
            }
        }
    }
    return; 
}
*/

/*
double DiagCal_det(const det& start)
{
    array<unsigned char,orbitalNum> occup = {0};
    for(int i=0; i< electronNum; i++)
        occup[start[i]] = 1;

    double sum = H.constant;

    int len = H.O_diag.size(); 
    for(int i=0; i<len; i++)
    {
        if(occup[H.O_diag[i][1]]!=0 ) sum += H.O_coeff_diag[i]; 
    } 

    len = H.OO_diag.size();
    for(int i=0; i<len; i++)
    {   
        int m = H.OO_diag[i][0];
        int n = H.OO_diag[i][1];
        int p = H.OO_diag[i][2];
        int q = H.OO_diag[i][3];

        if( (m ==n && p==q && occup[n]!=0 && occup[q]!=0) )
        {
            sum += 0.5*H.OO_coeff_diag[i];
        }   

        else if( (m ==q && n==p && occup[n]!=0 && occup[m]!=0))
        {
            sum -= 0.5*H.OO_coeff_diag[i];
        }
    }
    return sum;

}
*/

void ReadElements(Hamiltonian& H, string temp)
{
    string filename = temp;
    const char* u = &filename[0];
    FILE *fp = fopen(u,"r");

    //while(1)
    H.constant =0.0;

    for(int j = 0; j< 7; j++)
	{
        int num_line;
        fscanf(fp,"%d",&num_line);

        double coeff;
        int m,n,p,q;

        for(int i=0; i<num_line; i++) 
        {
            fscanf(fp, "%lf %d %d %d %d", &coeff, &m, &n, &p, &q);

            if(p==0 && q==0) 
            {
                if (m==n)
                {
                    if (m ==0 && n==0)
                    {
                        H.constant = coeff;
                    }
                    else
                    {   
                        oneBody temp;
                        temp[0] = m-2;
                        temp[1] = n-2;

                        H.O_diag.push_back(temp);
                        H.O_coeff_diag.push_back(coeff);
                    }
                }
                else 
                {
                    oneBody temp;
                    temp[0] = m-2;
                    temp[1] = n-2; 

                    H.O[temp[1]].push_back(temp);
                    H.O_coeff[temp[1]].push_back(coeff);
                }
            }
            else
            {
                twoBody temp;
                temp[0] = m-2;
                temp[1] = n-2;
                temp[2] = p-2;
                temp[3] = q-2;

                if ( (temp[0] == temp[1] && temp[2] == temp[3]) || (temp[0] == temp[3] && temp[2] == temp[1])  )
                {
                    ///diagonal element
                    H.OO_diag.push_back(temp);
                    H.OO_coeff_diag.push_back(coeff);
                }
                else
                {
                    if (!(m==n || n==p || p==q || m == q))
                    {
                        if (temp[1] < temp[3])
                        {
                            int idx = temp[1]+orbitalNum*temp[3];
                            H.OO[idx].push_back(temp);
                            H.OO_coeff[idx].push_back(coeff);
                        }
                    }
                    else
                    {
                        int index  = temp[1] < temp[3] ? temp[1]+orbitalNum*temp[3] : temp[3]+orbitalNum*temp[1];
                    /// locate the correct array

                        H.OO[index].push_back(temp);
                        H.OO_coeff[index].push_back(coeff);
                    }

                    /*
                    /// off-diagonal
                    int index  = temp[1] < temp[3] ? temp[1]+orbitalNum*temp[3] : temp[3]+orbitalNum*temp[1];
                    /// locate the correct array

                    H.OO[index].push_back(temp);
                    H.OO_coeff[index].push_back(coeff);
                    */
                }
            }
        }

        if(num_line == 1) break;
    }
    fclose(fp);
}


void load_integrals(string filename)
{
    cout.precision(12);
    /// initialize H
    H.OO.resize(orbitalNum*orbitalNum);
    H.OO_coeff.resize(orbitalNum*orbitalNum);

    H.O.resize(orbitalNum);
    H.O_coeff.resize(orbitalNum);

    ReadElements(H,filename);
    /*
    det start_det = FindHF();

	//{1,2,9,10,14,16,18,20,22,24,26};//{2,3,10,11,14,16,18,20,22,24,26};

    psi.basis.push_back(start_det);
    psi.amp.push_back(1.0);
    psi.pos[start_det] = 0;
    */
 }

/*
void printDet(const det& u)
{
    for(int i=0; i<electronNum; i++)
        printf("%hhu ",u[i]); 
    cout <<endl;
}
*/

unsigned int mergeBuffer(dynVec<DET_TYPE,FLOAT_TYPE>& buff, unsigned int buffSize)
{
    /*
        Must be sorted before merging,
    */

    if (buffSize <=1)
        return buffSize;  
        
    unsigned int mergedSize = 0; 
    DET_TYPE ptrBasis = buff.basis[0]; 
    FLOAT_TYPE ptrAmp = buff.amp[0];   

    for (int i=1; i<buffSize; i++)
    {
        if (ptrBasis != buff.basis[i])
        {
            buff.basis[mergedSize] = ptrBasis; 
            buff.amp[mergedSize] = ptrAmp; 
            
            ptrBasis = buff.basis[i];
            ptrAmp = buff.amp[i]; 

            mergedSize++; 
        }
        else 
        {
            ptrAmp += buff.amp[i]; 
        }
    }  
    // Watch out the tail 
    buff.basis[mergedSize] = ptrBasis; 
    buff.amp[mergedSize] = ptrAmp; 

    return mergedSize+1; 
}

void SpMVMtply_bucket(dynVec<DET_TYPE, FLOAT_TYPE>& wf_input_o, 
                      dynVec<DET_TYPE, FLOAT_TYPE>& wf_output, 
                      double lambda, double coeff, FLOAT_TYPE& diag_SpMSpV)
{   
    /*
        assume all elements have acsending order 
        How to allocate memory dynamically? 
        Calculate <v|H|v> here to avoid loss of accuracy
    */

    using colvec = dynVec<DET_TYPE,FLOAT_TYPE>; 

	if (wf_input_o.len ==0)
	{
		wf_output.setZeroVec();
		return;
	}

    dynVec<DET_TYPE, FLOAT_TYPE> wf_input(wf_input_o.cap);
    wf_input.naiveCopy(wf_input_o); 

    introsort_amp(wf_input,wf_input.len);
	
    unsigned int len = wf_input.len; 
    unsigned int expand_fac = bucket_unit;     
	unsigned int expand_fac_final = 5000; 
	
	unsigned int maxFinalLen = 10000000;

	// when buff is small no need to apply cutoff 
    unsigned int maxLocalBuff = maxFinalLen / (bucketNum*nThread)*2;
	//cout << maxLocalBuff <<endl; 
    unsigned int st_idx = wf_input.len < Nd ? 0 : wf_input.len - Nd;  
    unsigned int end_idx = wf_input.len; 

	unsigned int total_merged = 0;

	FLOAT_TYPE cutoff_ratio = 0.80; 
    FLOAT_TYPE cutoff_val = SpMSpVFinalResCut;

    FLOAT_TYPE diagAlpha = 0.0; 
    
    //Allocate space for buckets 
	vector<dynVec<DET_TYPE,FLOAT_TYPE>> resGlobal; 
    dynVec<DET_TYPE, FLOAT_TYPE> resFinal;  

    vector<colvec> wf_input_bucketed(bucketNum);

    resGlobal.resize(nThread*bucketNum);
    resFinal.reserveMem(len*expand_fac_final);

    /*
        idx = thread_idx + thread_num * bucket_idx
    */  
	
    vector<unsigned int> bucketCnt(bucketNum,0);

    for(int i = st_idx; i< end_idx; i++) 
    {
        auto bcktN = indexCal(wf_input.basis[i]) / bucketSize; 
        bucketCnt[bcktN] ++;         
	}

    /*
        initialize wf_input_bucketed
    */ 
	
    for (int bucketIdx =0; bucketIdx < bucketNum; bucketIdx ++)
    {
        //determinant number in each bucket 
        wf_input_bucketed[bucketIdx].reserveMem(bucketCnt[bucketIdx]);
    }

    for(int i = st_idx; i< end_idx; i++) 
    {
        //fill wf_input_bucketed   
        auto bcktN = indexCal(wf_input.basis[i]) / bucketSize; 
        wf_input_bucketed[bcktN].cpyFromPair(wf_input.basis[i], wf_input.amp[i]);
    }

    for (int bucketIdx =0; bucketIdx < bucketNum; bucketIdx ++)
    {
        //    make each bucket sorted 
        introsort(wf_input_bucketed[bucketIdx], wf_input_bucketed[bucketIdx].len);
    }
	


    for (int i = 0; i< nThread*bucketNum; i++)
    { 
		///no prediction
		
		//resGlobal[i].reserveMem((len/(nThread*bucketNum)+1) * expand_fac);      

		///with prediction
		int bckt_idx = i / nThread;
        resGlobal[i].reserveMem(((bucketCnt[bckt_idx] / nThread) + 4) * expand_fac);   
	}   

	cout << "bucket mem done\n";
	cout << st_idx  <<' '<<end_idx <<endl;
	auto t1 = system_clock::now();    
	 
    #pragma omp parallel    
    {          
        #pragma omp for schedule(static, 25) 
        for (int i=st_idx; i< end_idx; i++) 
        {   
            unsigned int thread_idx = omp_get_thread_num(); 
			
			auto basis_i = wf_input.basis[i];
            auto amp_i = wf_input.amp[i]; 

			double temp = coeff * amp_i;
			OffDiagGen_bit(basis_i,temp,resGlobal,thread_idx);  
            
            double diag_coeff = (DiagCal_bit(basis_i) - lambda) * temp;  
			unsigned int pos = thread_idx + nThread * (indexCal(basis_i) >> bucketSize_bit ); 
			resGlobal[pos].cpyFromPair(basis_i,diag_coeff);
        }

        #pragma omp barrier 

        /*  
		   SPA is much slower than sort-and-merge 
        */      
		
        #pragma omp for schedule(dynamic)
        for (int i=0; i<bucketNum; i++)
        {
            unsigned int localSize = 0;
            unsigned int mergedSize; 

            for (int j = 0; j < nThread; j++)
            {    
                unsigned int pos = j + nThread * i;  
                localSize += resGlobal[pos].len;
            }   
            
            dynVec<DET_TYPE, FLOAT_TYPE> localBuff;
            localBuff.reserveMem(localSize);        

            for (int j = 0; j < nThread; j++)
            {
                unsigned int pos = j + nThread * i;  
                localBuff.cpyFromVec(resGlobal[pos]);
            } 

            introsort(localBuff, localSize); 
			mergedSize = mergeBuffer(localBuff,localSize); 

            localBuff.len = mergedSize; 

            /*
                inner product in the ith bucket 
            */
            
			FLOAT_TYPE diag_temp = localBuff.dot(wf_input_bucketed[i]);

            #pragma omp atomic
            diagAlpha += diag_temp; 

			introsort_amp(localBuff,mergedSize); 

            #pragma omp critical
            {   
				/*
                    Copying the important determinants from localBuff     
					
					Two strategies:
						1. keep those with amplitude larger than a certain threashold 
						2. keep top 5% of each bucket
				*/  	



				total_merged += mergedSize; 
				//unsigned int localBuff_start = mergedSize *(1 - cutoff_ratio) < maxLocalBuff ? 0 : mergedSize*cutoff_ratio;
                //unsigned int localBuff_start = mergedSize < maxLocalBuff ? 0 : mergedSize - maxLocalBuff; 
				//unsigned int localBuff_end = mergedSize; 

				//unsigned int offset = resFinal.len - cutoff; // ?  			
                unsigned int copy_idx=0;	
                
                // strategy 1 
                for(unsigned int k = 0; k < mergedSize; ++k) 
				{
                    // taking in non-zero dets 
                    if (abs(localBuff.amp[k])>cutoff_val)
                    {    
                        resFinal.basis[resFinal.len + copy_idx] = localBuff.basis[k];
                        resFinal.amp[resFinal.len + copy_idx] = localBuff.amp[k];    
                        copy_idx++;
                    }
                    /*
                    else 
                    {
                        cout << localBuff.amp[k] << ' ' << abs(localBuff.amp[k]) <<endl;
                    }*/
                }
                

				/*
                // strategy 2
                for(unsigned int k = localBuff_end*cutoff_ratio; k < localBuff_end; ++k) 
                {
                    if (abs(localBuff.amp[k])>1e-14)
                    {    
                        resFinal.basis[resFinal.len + copy_idx] = localBuff.basis[k];
                        resFinal.amp[resFinal.len + copy_idx] = localBuff.amp[k];    
                        copy_idx++;
                    }
                }
				*/
                resFinal.len += copy_idx; 
                //(mergedSize - cutoff); 
            }   
		}
    }   
	cout<< "total merged:" << total_merged<<endl;
    cout <<"final raw vec len:"<< resFinal.len <<endl;
    cout << "value of exact alpha:" << diagAlpha<<endl;  
    
    diag_SpMSpV = diagAlpha; 

    /// stat for memory usage
	
	double r = 0; 
	for (int i=0; i< resGlobal.size(); i++)
    {
		//if (resGlobal[i].len > resGlobal[i].cap) 
        //cout << resGlobal[i].len << ' ' <<  resGlobal[i].cap <<endl; 
		double a = resGlobal[i].len;
		double b = resGlobal[i].cap;
		//cout << a <<' '<<b<<endl;
		a /= b;
		
		//cout << "occup rate:" << a <<endl<<endl;
		if (a>r) r = a; 
	}
	cout <<"max bucket occup ratio:" << r <<endl;
    
    
    auto t2 = system_clock::now();      

    
    unsigned int st = resFinal.len < wf_output.cap ? 0 : resFinal.len - wf_output.cap;
	
	//the following step removes dets with 0 amplitude
	int copy_idx =0 ;

    for(int i=st; i < resFinal.len; i++)
    {
        //if( abs (resFinal.amp[i]) > 1e-14)
		//{
		wf_output.basis[i-st] = resFinal.basis[i];  
		wf_output.amp[i-st] = resFinal.amp[i]; 
            //copy_idx ++;
		//}
		 
	}   
    wf_output.len = resFinal.len - st; 
	//wf_output.len = resFinal.len - st; 

    // put back to determinant-major 
	introsort(wf_output,wf_output.len);    

    auto t3 = system_clock::now();

    std::chrono::duration<double> d1 = t2-t1;
    std::chrono::duration<double> d2 = t3-t2;

    cout << "generating time: "<< d1.count() <<endl;
    cout << "merging time " << d2.count() <<endl; 
    
    return;
}
    
FLOAT_TYPE calEnergy(dynVec<DET_TYPE,FLOAT_TYPE>& psi)
{
    /*
        expectation energy  
    */  

    FLOAT_TYPE DiagCal_bit(const det_bit);
    void OffDiagGen_bit_e(const DET_TYPE, dynVec<DET_TYPE,FLOAT_TYPE>&); 

    FLOAT_TYPE energy = 0.0; 

    unordered_map<DET_TYPE,unsigned int> pos;
    pos.reserve(psi.len);

    for(unsigned int i=0; i<psi.len; i++)
        pos[psi.basis[i]] = i; 

    #pragma omp parallel for reduction(+ : energy)
    for (unsigned int i=0; i< psi.len; i++)
    {    
        auto current_det = psi.basis[i];
        auto current_amp = psi.amp[i];

        FLOAT_TYPE energy_i = current_amp * DiagCal_bit(current_det); 

        dynVec<DET_TYPE,FLOAT_TYPE> tgt; 
        OffDiagGen_bit_e(current_det,tgt);
        /// ? something is missing here    
        //cout << tgt.len <<endl;
		for (int j=0; j< tgt.len; j++)
        {
            if(pos.find(tgt.basis[j]) != pos.end())
            {
                energy_i += tgt.amp[j] * psi.amp[pos[tgt.basis[j]]]; 
            }
        }   
        energy += current_amp * energy_i;   
        
    }       
    return energy;  

}
