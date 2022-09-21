#include "all.h"

int spinCal_bit(const det_bit k)
{
	det_bit one = 1;

    int spin = 0;

    // occupied orbitals in increasing order
	for (int i=0;i<orbitalNum; i++)
	{
        if ( (k & (one << i)) !=0) 
        {
        	if (i%2 ==0 )  spin ++;
            else spin--; 
        } 
	}
 
	return spin; 
}	

/*
det bitToDet(const det_bit k)
{
	det u; 
	det_bit one = 1;

    int pos = 0;

    // occupied orbitals in increasing order
	for (int i=0;i<orbitalNum; i++)
	{
        if ( (k & (one << i)) !=0) 
        {
        	u[pos] = i; 
        	pos++; 
        } 
	}
	return u; 
}	
*/
/*
det_bit detToBit(const det k)
{
	det_bit result = 0; 
	det_bit one = 1; 
	for (int i=0; i< electronNum; i++)
	{
		result += (one << k[i]); 
	}
	return result;
}
*/
unsigned int nck(int n, int k)
{
    // k=0 or k=n -> 1 
    // 
    int kk = k > (n/2) ? (n-k) : k;
    int ans = 1;  

    for (unsigned int i =0 ; i< kk; i++)
    {
        ans *= (n-i);
        ans /= (i+1); 

    }
    return ans; 
}   

unsigned long long indexCal(const det_bit& ref)
{
    /// good for Sz =0; but can be easily extended to arbitrary Sz 
    /// i%2 ==0: spin up else spin down 
    
    det_bit one = 1;
    det_bit maskSpinUp = 0x5555555555555555;  // 5 : 0101
    //det_bit maskSpinDown = 0xAAAAAAAAAAAAAAAA; // A : 1010

    unsigned int idxUp = 0;
    //unsigned int idxDn = 0; 
    unsigned int posUp = 0; 
    //unsigned int posDn = 0;
             
    // remove the spin down bits
    maskSpinUp &= ref;  

    // last 1-bit trick 
    while(maskSpinUp) 
    {
        //get the orbital number with ctz 
        int orb = __builtin_ctzll(maskSpinUp) >> 1;

        idxUp += nck_list[posUp][orb]; 

        maskSpinUp = maskSpinUp & (maskSpinUp-1);
        
        posUp++;
    }   
    //__builtin_ctzll 
    // this loop is a waste
    /*
    for (unsigned int i=0; i< orbitalNum/2; i++)
    {
        unsigned long long tempUp = ref &(one << (2*i));
        //unsigned long long tempDown = ref &(one << (2*i+1));

        if (tempUp != 0 ) 
        {
            if (i != posUp) // means this electron has hopped 
            {   
                idxUp += nck_list[posUp][i]; 
            }  
            posUp++; 
        }
        
        //if (tempDown != 0)
        //{
        //    if (i != posDn) 
        //    {        
        //        idxDn += nck(i,posDn+1);//nck_list[posDn][i];//; 
        //    }
        //    posDn++;
        //}
    }
    */
    return idxUp;// + idxDn * blockDim;
}   
    
double signCal_bit(const det_bit& ref, int p, int q)
{
	//assume p < q. using bit tricks
    
	bool cnt = true; 
	det_bit temp = ref; 
	int one = 1; 

	temp >>= (p+1); 
	temp <<= (p+1);
	temp <<= (sizeof(ref)*8-q);  

    // count how many 1s are in the section 
    int t = __builtin_popcountll(temp); 

    return (t & one) == 1 ? -1.0 : 1.0;  
    // replace with builtin functions __builtin_popcount() 
	/*
    while(temp)	
	{
		temp = temp & (temp-1);
		cnt = !cnt; 
	}	
	return cnt == true ? 1.0 : -1.0; 
    */

    //__builtin_popcountll(temp) == 0 ? 
}	

void OffDiagGen_bit(const DET_TYPE reff, 
                    const FLOAT_TYPE coeff, 
					vector<dynVec<DET_TYPE,FLOAT_TYPE>>& tgt,
                    unsigned int thread_idx) 
{
    //void printDet(const det&);  

    det_bit one = 1;
    //array<unsigned char,electronNum> occup_orbs={0};
    vector<unsigned char> occup_orbs(electronNum, 0); 
    int pos = 0;    

    // occupied orbitals in increasing order
	for (int i=0;i<orbitalNum; i++)
	{
        if ( (reff & (one << i)) !=0) 
        {
        	occup_orbs[pos] = i; 
        	pos++; 
        }
	}
		
    /// two body operators
    for (int i=0; i<electronNum; i++)
        for(int j=i+1;j<electronNum; j++)
        {
            //unsigned idx = reff[i]  + reff[j] * orbitalNum; 
            // idx determines which two orbitals are occupied 
        	unsigned idx = occup_orbs[i] + occup_orbs[j] * orbitalNum; 
            unsigned len = H.OO[idx].size();
            
            for(int k = 0; k< len; k++)
            {
                det_bit start = reff; 

                unsigned char m = H.OO[idx][k][0];
                unsigned char n = H.OO[idx][k][1];
                unsigned char p = H.OO[idx][k][2];
                unsigned char q = H.OO[idx][k][3];
                
                ///double excitation	

               	auto m_stat = reff & (one << m);
               	auto p_stat = reff & (one << p);

                if ( (m!=n) && (n!=p) && (p!=q) && (m!= q) && (m_stat == 0 ) && (p_stat == 0) )
                {	
                    
                	double sign; 

                	sign = p < q ? signCal_bit(start, p, q) : signCal_bit(start, q, p);
                	start = start - (one<<q) + (one << p); 

                	sign =  sign * (m < n ? signCal_bit(start, m, n) : signCal_bit(start, n, m));
                	start = start - (one<<n) + (one << m); 

                    //unsigned long long start_idx = thread_idx + nThread * (indexCal(start) / bucketSize); 
                    unsigned long long start_idx = thread_idx + nThread * (indexCal(start) >> bucketSize_bit); 

                    //if(start_idx >= tgt.size())
                    //    cout << thread_idx <<' ' <<nThread << ' ' << (indexCal(start) >> bucketSize_bit)<<endl; 
					tgt[start_idx].cpyFromPair(start,coeff*sign*H.OO_coeff[idx][k]); 
                    
                    //tgt.cpyFromPair(start,coeff*sign*H.OO_coeff[idx][k]); 

                }
                                
                ///single excitation
                else if( m!=n && n!=p && p==q && q!=m && m_stat == 0 )
                {
                    double sign; 

                	sign =  m < n ? signCal_bit(start, m, n) : signCal_bit(start, n, m);
                	start = start - (one<<n) + (one << m); 
                
                    //unsigned long long start_idx = thread_idx + nThread*(indexCal(start) / bucketSize);
                    unsigned long long start_idx = thread_idx + nThread*(indexCal(start) >> bucketSize_bit ); 

                    //if(start_idx >= tgt.size())
                    //    cout << thread_idx <<' ' <<nThread << ' ' << (indexCal(start) >> bucketSize_bit)<<endl; 

					tgt[start_idx].cpyFromPair(start,coeff*sign*H.OO_coeff[idx][k]); 
                    
                    //tgt.cpyFromPair(start,coeff*sign*H.OO_coeff[idx][k]); 
                }

                else if ( m!=n && n==p && p!=q && q!=m && m_stat == 0 )
                {
                	double sign; 
                    
                	sign =  q < m ? signCal_bit(start, q, m) : signCal_bit(start, m, q);
                	start = start - (one<<q) + (one << m); 
				   
                    //unsigned long long start_idx = thread_idx + nThread *(indexCal(start) / bucketSize); 
                    unsigned long long start_idx = thread_idx + nThread *(indexCal(start) >> bucketSize_bit); 
                    //if(start_idx >= tgt.size())
                    //    cout << thread_idx <<' ' <<nThread << ' ' << (indexCal(start) >> bucketSize_bit)<<endl; 
					tgt[start_idx].cpyFromPair(start,-coeff*sign*H.OO_coeff[idx][k]); 
                    
                    //tgt.cpyFromPair(start,coeff*sign*H.OO_coeff[idx][k]);
                }

            }

        }		
    /// one body operators

    for (int i=0; i<electronNum; i++)
    {
        unsigned idx = occup_orbs[i];
        unsigned len = H.O[idx].size();

        for (int j = 0; j< len; j++)
        {
            det_bit start = reff; 
            unsigned char m = H.O[idx][j][0];
            unsigned char n = H.O[idx][j][1];
            
            auto m_stat = reff & (one << m); 
            
            if (m_stat == 0 ) 
            {
                double sign; 
                sign =  m < n ? signCal_bit(start, m, n) : signCal_bit(start, n, m);
                start = start - (one<<n) + (one << m); 
                
                //unsigned long long start_idx = thread_idx + nThread*(indexCal(start) / bucketSize); 
                unsigned long long start_idx = thread_idx + nThread*(indexCal(start) >> bucketSize_bit ); 
				
				/*
                    if(start_idx >= tgt.size())
                        cout << thread_idx <<' ' <<nThread << ' ' << (indexCal(start) >> bucketSize_bit)<<endl; 
				*/
                tgt[start_idx].cpyFromPair(start,coeff*sign*H.O_coeff[idx][j]); 
                    
                //tgt.cpyFromPair(start,coeff*sign*H.O_coeff[idx][j]); 
            }
        }
    }   
    return; 
}

void OffDiagGen_bit_e(const DET_TYPE reff, dynVec<DET_TYPE,FLOAT_TYPE>& tgt)
                     
{
    tgt.reserveMem(numConnected); 

    det_bit one = 1;
    vector<unsigned char> occup_orbs(electronNum, 0); 
    int pos = 0;    

    // occupied orbitals in increasing order
    for (int i=0;i<orbitalNum; i++)
    {
        if ( (reff & (one << i)) !=0) 
        {
            occup_orbs[pos] = i; 
            pos++; 
        }
    }
    
    /// two body operators
    for (int i=0; i<electronNum; i++)
        for(int j=i+1;j<electronNum; j++)
        {
            //unsigned idx = reff[i]  + reff[j] * orbitalNum; 
            // idx determines which two orbitals are occupied 
            unsigned idx = occup_orbs[i] + occup_orbs[j] * orbitalNum; 
            unsigned len = H.OO[idx].size();
            
            for(int k = 0; k< len; k++)
            {
                det_bit start = reff; 

                unsigned char m = H.OO[idx][k][0];
                unsigned char n = H.OO[idx][k][1];
                unsigned char p = H.OO[idx][k][2];
                unsigned char q = H.OO[idx][k][3];
                
                ///double excitation    

                auto m_stat = reff & (one << m);
                auto p_stat = reff & (one << p);

                if ( (m!=n) && (n!=p) && (p!=q) && (m!= q) && (m_stat == 0 ) && (p_stat == 0) )
                {   
                    
                    double sign; 

                    sign = p < q ? signCal_bit(start, p, q) : signCal_bit(start, q, p);
                    start = start - (one<<q) + (one << p); 

                    sign =  sign * (m < n ? signCal_bit(start, m, n) : signCal_bit(start, n, m));
                    start = start - (one<<n) + (one << m); 

                    //unsigned long long start_idx = thread_idx + nThread * (indexCal(start) >> bucketSize_bit); 
                    //tgt[start_idx].cpyFromPair(start,coeff*sign*H.OO_coeff[idx][k]); 
                    
                    tgt.cpyFromPair(start,sign*H.OO_coeff[idx][k]); 

                }
                                
                ///single excitation
                else if( m!=n && n!=p && p==q && q!=m && m_stat == 0 )
                {
                    double sign; 

                    sign =  m < n ? signCal_bit(start, m, n) : signCal_bit(start, n, m);
                    start = start - (one<<n) + (one << m); 
                
                    //unsigned long long start_idx = thread_idx + nThread*(indexCal(start) / bucketSize);
                    //unsigned long long start_idx = thread_idx + nThread*(indexCal(start) >> bucketSize_bit ); 
                    //tgt[start_idx].cpyFromPair(start,coeff*sign*H.OO_coeff[idx][k]); 
                    
                    tgt.cpyFromPair(start,sign*H.OO_coeff[idx][k]); 
                }

                else if ( m!=n && n==p && p!=q && q!=m && m_stat == 0 )
                {
                    double sign; 
                    
                    sign =  q < m ? signCal_bit(start, q, m) : signCal_bit(start, m, q);
                    start = start - (one<<q) + (one << m); 
                   
                    //unsigned long long start_idx = thread_idx + nThread *(indexCal(start) / bucketSize); 
                    //unsigned long long start_idx = thread_idx + nThread *(indexCal(start) >> bucketSize_bit); 
                    //tgt[start_idx].cpyFromPair(start,coeff*sign*H.OO_coeff[idx][k]); 
                    
                    tgt.cpyFromPair(start,-sign*H.OO_coeff[idx][k]);
                }

            }

        }
    
    /// one body operators

    for (int i=0; i<electronNum; i++)
    {
        unsigned idx = occup_orbs[i];
        unsigned len = H.O[idx].size();

        for (int j = 0; j< len; j++)
        {
            det_bit start = reff; 
            unsigned char m = H.O[idx][j][0];
            unsigned char n = H.O[idx][j][1];
            
            auto m_stat = reff & (one << m); 
            
            if (m_stat == 0 ) 
            {
                double sign; 
                sign =  m < n ? signCal_bit(start, m, n) : signCal_bit(start, n, m);
                start = start - (one<<n) + (one << m); 
                
                //unsigned long long start_idx = thread_idx + nThread*(indexCal(start) / bucketSize); 
                //unsigned long long start_idx = thread_idx + nThread*(indexCal(start) >> bucketSize_bit ); 
                //tgt[start_idx].cpyFromPair(start,coeff*sign*H.O_coeff[idx][j]); 
                    
                tgt.cpyFromPair(start,sign*H.O_coeff[idx][j]); 
            }
        }
    }   
    return; 
}

FLOAT_TYPE DiagCal_bit(const det_bit reff)
{
    det_bit one = 1;
	double sum = H.constant;
    
    int len = H.O_diag.size(); 
    for(int i=0; i<len; i++)
    {
        if(( reff & (one << H.O_diag[i][1])) !=0 ) sum += H.O_coeff_diag[i]; 
    } 

    len = H.OO_diag.size();

    for(int i=0; i<len; i++) 
    {   
        unsigned m = H.OO_diag[i][0];
        unsigned n = H.OO_diag[i][1];
        unsigned p = H.OO_diag[i][2];
        unsigned q = H.OO_diag[i][3];

        /// 0 ~ unoccupied, otherwise occupied 
        auto isOccup_n = reff & (one << n);
        auto isOccup_q = reff & (one << q); 

        if (isOccup_n !=0 && isOccup_q !=0)
        {
        	if (m ==n && p==q) sum += 0.5*H.OO_coeff_diag[i];
        	else if (m ==q && n==p) sum -= 0.5*H.OO_coeff_diag[i];
        }

        /*
        if( (m ==n && p==q && isOccup_n !=0 && isOccup_q !=0) )
        {
            
        }

        else if( (m ==q && n==p && isOccup_n !=0 && isOccup_q!=0 )
        {
            
        }*/
    }
    return sum;

}

