 #include "all.h"

using colvec = dynVec<DET_TYPE,FLOAT_TYPE>; 

FLOAT_TYPE calEnergy(dynVec<DET_TYPE,FLOAT_TYPE>&);

//int spinNumCal(const det);
//det bitToDet(const det_bit); 

void calSD(DET_TYPE hf_det, colvec& result) 
{
    /**
        Diagonalize the single and double excitation block. 
        Similar to Semi-stochastic qmc and i-FCIQMC
    */
    void OffDiagGen_bit_e(const DET_TYPE, dynVec<DET_TYPE,FLOAT_TYPE>&);
    FLOAT_TYPE DiagCal_bit(const DET_TYPE); 
    void introsort(dynVec<DET_TYPE,FLOAT_TYPE>&, int); 

    colvec offDiag;
    OffDiagGen_bit_e(hf_det,offDiag);
    //cout << offDiag.len<<endl;
    vector<DET_TYPE> basis; 

    unordered_set<DET_TYPE> offDiag_set;
    unordered_map<DET_TYPE,int> pos; 

    basis.push_back(hf_det);
	
    for(int i=0; i<offDiag.len; i++)
    {
        if(offDiag_set.find(offDiag.basis[i]) == offDiag_set.end())
        {
            offDiag_set.insert(offDiag.basis[i]);
            basis.push_back(offDiag.basis[i]); 
        }  
    }
	cout << offDiag_set.size() <<endl; 
    int len = basis.size(); 
    cout << "SD block size: " << basis.size()<<endl; 

    for(int i=0; i<basis.size(); i++)
        pos[basis[i]] = i; 

    // A dense H 
    M HH(len,len);
    HH.setZero();

    /// fill the matrix 
    for(int i=0; i< len; ++i)
    {
        colvec pool; 
        OffDiagGen_bit_e(basis[i],pool);

        unordered_map<DET_TYPE,FLOAT_TYPE> pool_uom; 
        for (int j=0; j<pool.len; j++)
        {
            if (pool_uom.find(pool.basis[j]) ==  pool_uom.end() )
                pool_uom[pool.basis[j]] = pool.amp[j];
            else 
                pool_uom[pool.basis[j]] += pool.amp[j];
        }

        HH(i,i) = DiagCal_bit(basis[i]);

        for(auto ite = pool_uom.begin(); ite != pool_uom.end(); ite++)
        {
            auto det_now = ite->first;
            if(pos.find(det_now) != pos.end())
            {
                int p = pos[det_now];
                HH(i,p) = ite->second; 
            }
        }
        /*  
        for(int j=0; j<pool.len; j++)
        {
            /// only include those within the block
            if(pos.find(basis[j])!= pos.end())
            {
                int p = pos[pool.basis[j]];
                HH(i,p) += pool.amp[j];   
            }
        }*/
    }           


    SAES saes; 
    saes.compute(HH);
    auto eigenval = saes.eigenvalues();
    auto eigenvec = saes.eigenvectors();
    auto coeff = eigenvec.col(0);
	
	//cout << coeff <<endl; 

    cout <<"From CISD: "<<eigenval(0)<<endl<<endl;

    for(int i=0; i<len; i++)
    {
        result.cpyFromPair(basis[i],coeff[i]);
    }

	introsort(result,result.len);
}

void isOrdered(colvec& psi)
{
    for (int i=1; i<psi.len; i++)
    {
        if (psi.basis[i] <= psi.basis[i-1])
        {
            cout << "wrong\n";
            return; 
        }
    }
    cout << "Good\n";
}
/*
void checkPresence(colvec& src, colvec& tgt, FLOAT_TYPE alpha)
{

    unordered_map<DET_TYPE,FLOAT_TYPE> tgt_pool; 
    tgt_pool.reserve(tgt.len);  

    for(int i=0; i< tgt.len; i++)
        tgt_pool[tgt.basis[i]] = tgt.amp[i];

    FLOAT_TYPE max_amp = 0.0;

    for(int i=0; i< src.len; i++)
    {
        if (tgt_pool.find(src.basis[i]) != tgt_pool.end())
        {
            if (abs(alpha*tgt_pool[src.basis[i]] - src.amp[i])>1e-12)
            //cout << src.basis[i] <<endl;
                cout << alpha*tgt_pool[src.basis[i]] - src.amp[i]<<endl; 
        }   
    }

    cout << "max amplitude:" << max_amp <<endl; 
}

void isDuplicated(colvec& v)
{
    unordered_map<DET_TYPE,int> chi(v.len);

    for(int i=0; i< v.len; i++)
    {
        if (chi.find(v.basis[i]) != chi.end())
        {

            ;
        }
    }
}

*/
void spmspv_uom(colvec& v, double alpha)
{
    void OffDiagGen_bit_e(const DET_TYPE, dynVec<DET_TYPE,FLOAT_TYPE>&);
    FLOAT_TYPE DiagCal_bit(const DET_TYPE); 

    unordered_map<DET_TYPE,double> output(v.len*500);

    for(int i=0;i<v.len; i++)
    {
        colvec temp; 
        
        OffDiagGen_bit_e(v.basis[i],temp);
        temp.cpyFromPair(v.basis[i],DiagCal_bit(v.basis[i]));

        temp.scalarMtply(v.amp[i]);

        for(int j=0; j<temp.len; j++)
        {
            if (output.find(temp.basis[j]) == output.end())
                output[temp.basis[j]] = temp.amp[j];
            else 
                output[temp.basis[j]] += temp.amp[j];
        }
    }
    // output = H |psi_0> 
    cout <<"uom size:" << output.size() <<endl; 
    /*
    for(int i=0; i< tgt.len; i++)
    {
        if (abs(tgt.amp[i] -output[tgt.basis[i]]) >     1e-12)  
            cout <<  tgt.amp[i] <<' ' << output[tgt.basis[i]]<<endl; 
    }
    */

    for (int i=0; i<v.len; i++)
    {
        if (abs(alpha*v.amp[i] - output[v.basis[i]] ) >1e-12) 
            cout << alpha*v.amp[i] <<' '<< output[v.basis[i]]<<endl; 
    }
}
/*
void getBasis(colvec& v)
{
    int len_large = 20000000; 
    
    int numIterations = 5;

    vector<colvec> KrylovBasis;
    KrylovBasis.resize(numIterations);
    for(int i=0; i< numIterations; i++)
        KrylovBasis[i].reserveMem(Nd);

    KrylovBasis[0].naiveCopy(v);

    for (int i=1; i< numIterations; i++)
    {
        colvec u(len_large);
        SpMVMtply_bucket(KrylovBasis[i-1],u,0.0,1.0);
        KrylovBasis[i].shrinkFrom(u);
        KrylovBasis[i].normalize();
    }

    M S(numIterations,numIterations);
    S.setZero(); 

    for(int i=0; i< numIterations; i++)
        {

            for(int j=i; j< numIterations; j++)
            {
                S(i,j) = KrylovBasis[i].dot(KrylovBasis[j]);
                S(j,i) = S(i,j);
            }  
        }

    cout << S <<endl; 
}
*/
void calGroundEnergy(colvec& v,int numI,int numR)
{
    /**
        Core of Lanzcos.
        Routine A(2,7) by C. C. Paige.
    */
    void introsort(colvec&, int);
    void introsort_amp(colvec&, int);
    unsigned int mergeBuffer(colvec&, unsigned int);

    FLOAT_TYPE calEnergy(dynVec<DET_TYPE,FLOAT_TYPE>&);

    //cout << "num of threads used:" << nbThreads() <<endl<<endl;

    int numIterations = numI;
    int numRestarts = numR;

    int len_large = HPsiCap; 
    FLOAT_TYPE gamma;

    FLOAT_TYPE deltaE = 0.00001;
    FLOAT_TYPE previousE = 0.0; 
    FLOAT_TYPE nearestE = -1000000.0;
    
    int nearestIdx = 0; 

    //M KrylovBasis(spaceDim,numIterations);
    vector<colvec> KrylovBasis;
    KrylovBasis.resize(numIterations);
    for(int i=0; i< numIterations; i++)
        KrylovBasis[i].reserveMem(Nd);

    M KrylovM(numIterations,numIterations);
    KrylovM.setZero();
    
    M S(numIterations,numIterations);
    S.setZero(); 

    SAES saes;
    GSAES gsaes;

    colvec vv(numIterations*Nd); 

    for(int res = 0; res < numRestarts; res++)
    {   
		cout << "\n*****************block "; 
        cout << res;
        cout << " starts*****************\n";
	 
		KrylovBasis[0].naiveCopy(v); 

        colvec u(len_large);
        colvec w(Nd);
        colvec v_next(Nd);

        FLOAT_TYPE alpha_SpMSpV;

        SpMVMtply_bucket(v,u,0.0,1.0,alpha_SpMSpV);

        //u = H*|psi_0> 

		cout <<"H * phi_0 len:" << u.len <<endl; 

        //FLOAT_TYPE alpha = v.dot(u);
        FLOAT_TYPE alpha = alpha_SpMSpV; 
        //cout << "alpha with truncation:" << alpha <<endl; 

		KrylovM(0,0) = alpha; 

        FLOAT_TYPE beta;
        
		//spmspv_uom(v,alpha); 
        //checkPresence(u,v,alpha);

		/*
            <v|H|v> = alpha <v|v> does not necessarily mean |H|v> - alpha|v> = 0! 
            Use addTwoCorrected instead 
        */
        
        colvec temp(u.cap); 
        colvec temp_u(u.cap);

        cout << "*****************ROUND " << '0' << " FINISHED*****************\n";   

        for (int i=1; i<numIterations; i++)
        {   
            temp.naiveCopy(u);

            temp.addTwo(v,-alpha);
            //temp.addTwoCorrected(v,-alpha);
            
            v_next.shrinkFrom(temp);     

            cout << v_next.len <<endl;

            double gamma_trunc = v_next.norm(); 

            /// choose which beta to use 
            beta = gamma_trunc; 
            v_next.scalarMtply(1.0/gamma_trunc);
            
            cout << "matrix multiplication starts\n";

            //v_next contains the most important ones from full length |phi_i>  
            SpMVMtply_bucket(v_next,temp_u,0.0,1.0, alpha_SpMSpV);    
            
			cout << "matrix multiplication done\n";

            //alpha = v_next.dot(temp_u);
            alpha = alpha_SpMSpV; 
            KrylovM(i,i) = alpha;

            //cout << "alpha with truncation:" << alpha <<endl; 
            
            // calculate effective H matrix element
            for(int j= 0; j< i; j++)
            {
                KrylovM(j,i) = KrylovBasis[j].dot(temp_u);
                KrylovM(i,j) = KrylovM(j,i);
            }

            temp_u.addTwo(v,-beta);
            //temp_u.addTwoCorrected(v,-beta); 
            u.naiveCopy(temp_u);

            //alpha = v_next.dot(u);

            v.naiveCopy(v_next);
            KrylovBasis[i].naiveCopy(v);

            /*
                Filling the matrix
            */

            /*
            KrylovM(i,i-1) = KrylovBasis[i].dot(HPsi[i-1]);

            KrylovM(i-1,i) = KrylovM(i,i-1);
            KrylovM(i,i) = alpha;

            ///fill the whole effective H

            for(int j = 0 ; j< i-1; j++)
            {
                KrylovM(j,i) = KrylovBasis[i].dot(HPsi[j]); 
                KrylovM(i,j) = KrylovM(j,i); 
            }   
            */


            temp.setZeroVec();
            temp_u.setZeroVec();

            cout << "*****************ROUND " << i << " FINISHED*****************\n";

        }


        for(int i=0; i< numIterations; i++)
        {
            //cout << KrylovBasis[i].len <<endl;
            //for (int j = i; (j < i+2) && (j < numIterations); j++)
            for(int j=i; j< numIterations; j++)
            {
                S(i,j) = KrylovBasis[i].dot(KrylovBasis[j]);
                S(j,i) = S(i,j);
            }  
        }

        //cout << S <<endl<<endl; 
        
        gsaes.compute(KrylovM,S);
        auto lambda_gsaes = gsaes.eigenvalues();
        cout << "e from gsaes:\n" << lambda_gsaes <<endl;  

        auto eigenvec = gsaes.eigenvectors();

        /*
        
            filter out the spurious solutions 
        
        */
        
        double maxGap; 

        for(int i=0; i< lambda_gsaes.size(); i++)
        {
            if (i==0)
            {
                maxGap = abs(nearestE - lambda_gsaes(i));
                
				nearestIdx = i; 
            }
            else
            {
                if (abs(nearestE - lambda_gsaes(i)) < maxGap)
                {
                    maxGap = abs(nearestE - lambda_gsaes(i)); 
                    nearestIdx = i; 
                }
            }
        }
		nearestE = lambda_gsaes(nearestIdx);
        
        auto coeff = eigenvec.col(nearestIdx);

        cout << "coeff from gsaes:\n" << coeff << endl; 
        cout << "coeff norm:" << coeff.norm() <<endl; 		
        cout << "vBv:" << coeff.transpose() * S * coeff<<endl;
        cout << "vHv:" << coeff.transpose() * KrylovM * coeff <<endl; 

        /*
            Hv = eSv 
        */
        
        int len_vv = 0; 

        for(int i=0; i< numIterations; i++)
            len_vv += KrylovBasis[i].len;

		colvec vv(len_vv);

        for(int i=0; i< numIterations; i++)
        {
            KrylovBasis[i].scalarMtply(coeff(i));
            vv.cpyFromVec(KrylovBasis[i]);
        }   

        introsort(vv,vv.len);
        unsigned int merged_len = mergeBuffer(vv,vv.len);
		
        cout <<"merged len:"<< merged_len <<endl; 
        vv.len = merged_len;

        FLOAT_TYPE vv_norm = vv.norm();
        cout << "vv norm:" << vv_norm <<endl; 

        //cout << "energy vv before truncation:" << calEnergy(vv)/vv_norm<<endl;

        v.shrinkFrom(vv);
        v.normalize(); 
        
        cout << "energy vv after truncation:" << calEnergy(v)<<endl;

		cout << KrylovM <<endl; 		

        KrylovM.setZero();
        for (int i=0; i< numIterations; i++)
            KrylovBasis[i].setZeroVec();    
        
        //cout << "residual:" <<residual<<endl;
        cout << "gsaes energy:" << lambda_gsaes(nearestIdx) <<endl; 


        cout << "*****************block "; 
        cout << res;
        cout << " ends*****************\n";

        if (previousE - lambda_gsaes(nearestIdx) < deltaE)
        {
            break; 
        }
        previousE = lambda_gsaes(nearestIdx); 
        //if (residual < threshold )
        //    break;
    }
	cout << "final e:"<< calEnergy(v)<<endl;

}
