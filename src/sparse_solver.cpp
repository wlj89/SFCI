 #include "all.h"

/*
    sparse solvers
    
*/

using colvec = dynVec<DET_TYPE,FLOAT_TYPE>; 



//int spinNumCal(const det);
//det bitToDet(const det_bit); 

void calc_CISD(DET_TYPE hf_det, colvec& result, int state_idx) 
{
    /**
        Diagonalize the CISD block. 

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

    auto coeff = eigenvec.col(state_idx);
	
    cout << "10 lowest states" <<endl; 

    cout << eigenval.head(10) <<endl;
	//cout << coeff <<endl; 

    cout <<"From CISD: "<<eigenval(state_idx)<<endl<<endl;

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
void sparse_krylov_solve_eigenval(colvec& v,int numI,int numR)
{
    /**
        Core of Lanzcos.
        Routine A(2,7) by C. C. Paige.
    */

    using time_dur = std::chrono::duration<double>; 

    void introsort(colvec&, int);
    void introsort_amp(colvec&, int);
    unsigned int merge_bucket(colvec&, unsigned int);

    FLOAT_TYPE calc_energy(const dynVec<DET_TYPE,FLOAT_TYPE>&);

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
            auto t_res_start = system_clock::now();

            temp.naiveCopy(u);

            temp.addTwo(v,-alpha);
            //temp.addTwoCorrected(v,-alpha);
            

            /// why do we need the shrink from??? 
            v_next.shrinkFrom(temp);     

            cout << v_next.len <<endl;

            double gamma_trunc = v_next.norm(); 

            /// choose which beta to use 
            beta = gamma_trunc; 
            v_next.scalarMtply(1.0/gamma_trunc);
            
            cout << "matrix multiplication starts\n";

            auto t_spmspvm_start = system_clock::now();
            //v_next contains the most important ones from full length |phi_i>  
            SpMVMtply_bucket(v_next,temp_u,0.0,1.0, alpha_SpMSpV);    
            auto t_spmspvm_end = system_clock::now();

            time_dur d_spmspvm = t_spmspvm_end - t_spmspvm_start;

            cout << "Single spmspvm time:" << d_spmspvm.count() <<endl;
			cout << "matrix multiplication done\n";

            //alpha = v_next.dot(temp_u);
            alpha = alpha_SpMSpV; 
            KrylovM(i,i) = alpha;

            //cout << "alpha with truncation:" << alpha <<endl; 
            
            // calculate effective H matrix element
            //auto t_Hessenberg_start = system_clock::now();

            for(int j= 0; j< i; j++)
            {
                KrylovM(j,i) = KrylovBasis[j].dot(temp_u);
                KrylovM(i,j) = KrylovM(j,i);
            }

            //auto t_Hessenberg_end= system_clock::now();
            //time_dur d_Hess = t_Hessenberg_end-t_Hessenberg_start;
            //cout << "time of filling Hessenberg Matrix:" << d_Hess.count() <<endl;

            auto t_addTwo_start = system_clock::now();

            temp_u.addTwo(v,-beta);

            auto t_addTwo_end = system_clock::now();
            time_dur d_addTwo = t_addTwo_end - t_addTwo_start;
            cout << "time of addTwo function:" << d_addTwo.count()<<endl; 
            
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

            auto t_res_end = system_clock::now();
            time_dur d_res = t_res_end - t_res_start;
            cout << "time of current iter:" << d_res.count() << endl;

            cout << "*****************ROUND " << i << " FINISHED*****************\n";

        }

        //auto t_S_start = system_clock::now();

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
        //auto t_S_end = system_clock::now(); 

        //time_dur d_S = t_S_end - t_S_start;
        //cout << "filling the S matrix:" << d_S.count()<<endl; 

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
        unsigned int merged_len = merge_bucket(vv,vv.len);
		
        cout <<"Psi_res len:"<< merged_len <<endl; 
        vv.len = merged_len;

        FLOAT_TYPE vv_norm = vv.norm();
        cout << "Psi_res norm:" << vv_norm <<endl; 

        //cout << "energy vv before truncation:" << calEnergy(vv)/vv_norm<<endl;

        v.shrinkFrom(vv);
        v.normalize(); 
        
        cout << "energy Psi_res after truncation:" << calc_energy(v)<<endl;

		//cout << KrylovM <<endl; 		

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
	cout << "final e:"<< calc_energy(v)<<endl;


}


void sparse_krylov_solve_eigenval_2(colvec& v,int numI,int numR)
{
    /*
        a more efficient sparse eigen solver
    */

    
    
}

void sparse_krylov_solve_linear(colvec& b, 
                                colvec& result_final,
                             FLOAT_TYPE target_energy,
                           unsigned int num_iter,
                           unsigned int num_restart)
{
    /*
        sparse Krylov linear solver under the SFCI wavefunction.
        use x0 = 0  and b with |b|=1  
    */
    void spmspvm_bucket_linear( const colvec&,
                                      colvec&,
                                   FLOAT_TYPE, 
                          vector<FLOAT_TYPE>&,
                      vector<vector<colvec>>&,
                                      int,int);

    int make_bucket(const colvec&, vector<colvec>&);

    void add_bucketed_vec(vector<vector<colvec>>&,
                                          colvec&,
                                          ColVec&);

    FLOAT_TYPE spmspvm_bucket_res_norm( const colvec&,
                                             colvec&,
                                        vector<colvec>&,
                                            FLOAT_TYPE,
                                                int);

    FLOAT_TYPE DiagCal_bit(const DET_TYPE);

    int max_bkt_size=0; 

    //M hessenberg_rec(num_iter+1,num_iter);
    M hessenberg(num_iter+1, num_iter+1);
    M s(num_iter+1,num_iter+1);

    vector<FLOAT_TYPE> Hij(num_iter+1);

    //basis
    vector<colvec> krylov_basis;
    krylov_basis.resize(num_iter+1);

    for (int i=0; i< num_iter+1; i++)
        krylov_basis[i].reserveMem(Nd); 

    //buckted basis
    vector<vector<colvec>> krylov_basis_bucketed;
    krylov_basis_bucketed.resize(num_iter+1);

    vector<colvec> b_buckted;
    make_bucket(b,b_buckted);

    colvec result;  // Vm*y
    colvec w(Nd); 
    colvec r0(Nd);

    colvec x_m(Nd);

    // explicitly assume x0=0
    r0.naiveCopy(b);

    for (int res = 0; res < num_restart; res++)
    {
        double beta = r0.norm(); 

        krylov_basis[0].naiveCopy(r0);
        krylov_basis[0].normalize();
   
        cout << "initial residue: " << beta <<endl;  

        for (int i=0; i< num_iter+1; i++)   
        {
            int bkt_size_temp = make_bucket(krylov_basis[i], krylov_basis_bucketed[i]);
            
            //cout << "basis vector bucketed" <<endl; 

            max_bkt_size = max(max_bkt_size,bkt_size_temp);

            spmspvm_bucket_linear(  krylov_basis[i],
                                    w,
                                    target_energy,
                                    Hij,
                                    krylov_basis_bucketed,
                                    max_bkt_size,
                                    i);

            for(int j=i; j>=0; j--)
            {
                hessenberg(j,i) = Hij[j];
                hessenberg(i,j) = hessenberg(j,i);

                s(i,j) = krylov_basis[j].dot(krylov_basis[i]);
                s(j,i) = s(i,j);  
            }

            w.normalize(); 

            if ( i< num_iter )
                krylov_basis[i+1].naiveCopy(w);

        }   

        M hessenberg_rec = hessenberg.block(0,0,num_iter+1,num_iter);

        //cout << hessenberg_rec <<endl; 

        cout << "hessenberg\n"  << hessenberg <<endl <<endl; 

        cout << "overlap\n" << s <<endl; 

        /*
            solve the least sqaure 

            or, use inverse as in FOM? 
        */

        ColVec e1 = s.col(0);
        e1 = e1*beta;

        ColVec y = (s * hessenberg_rec).colPivHouseholderQr().solve(e1);

        /*
            calculate xm = x0+ Vm*y 
        */  
        add_bucketed_vec(krylov_basis_bucketed, result, y);

        cout << "length of final result: ";
        cout << result.len << endl;
        //colvec result_shrinked(Nd);

        /*
            calculate |b - Ax|
        */
        /*
        cout << "actual residue: " << spmspvm_bucket_res_norm( result, 
                                                        b_buckted, 
                                                        target_energy,
                                                        max_bkt_size) <<endl;
        
        */
        /*
            x_m = x_m + Vm *y 
        */
        x_m.addTwo(result,1.0);
        //x_m.shrinkFrom(result);

        // also return r0 = b - A*x_m 
        cout << "residue after truncation: " 
             << spmspvm_bucket_res_norm( x_m, 
                                         r0,
                                         b_buckted, 
                                         target_energy,
                                         max_bkt_size) <<endl;

        /*
            returned r0 is in fact A*x_m - b. 
            need a minus sign. 
        */
        r0.scalarMtply(-1.0);
        //cout <<endl;
        cout << "norm of r0:" << r0.norm() <<endl;

        cout << "*************************restart " << res << " ends*************************\n"; 
    }
    
    result_final.naiveCopy(x_m);

    return; 

}

void sparse_solve_SI(colvec& b, FLOAT_TYPE target_energy)
{

    /*

        lambda = (Ei + E_(i-1)) /2
    */

    FLOAT_TYPE calc_energy(const colvec&);

    colvec result(Nd);

    double Ei = target_energy;

    for (int i=0; i< 1; i++ )
    {
        cout << "################################################################\n";
        cout << i << "th SI starts"  <<endl; 
        cout << "using target energy: " << target_energy <<endl;

        sparse_krylov_solve_linear(b,result,target_energy,numIteration,numRestart); 

        b.naiveCopy(result);
        b.normalize();

        FLOAT_TYPE Ei_temp =  calc_energy(b);
        cout << "energy after " << i << "th SI:" << Ei_temp <<endl;


        if (abs(Ei_temp - Ei) < 5E-6)  break;

        Ei = Ei_temp; 

        //target_energy = 0.5*(target_energy + Ei);

        cout << "################################################################\n";
    }

    return;
}

void sparse_main(FLOAT_TYPE target_energy)
{
    void calc_CISD(DET_TYPE, colvec&, int);
    /*
    void generateBasis( vector<DET_TYPE>&, unordered_map<DET_TYPE,int>&);
    
    unsigned seed = 6123123;

    colvec sparse_b(Nd);

    vector<DET_TYPE> basis;    
    unordered_map<DET_TYPE,int> basis_idx;

    generateBasis(basis,basis_idx);
    int spaceDim = basis.size();

    cout << "space dim:" << spaceDim <<endl;
    sparse_b.cpyFromPair(basis[spaceDim/223+61],1.0);
    */

    colvec b0(Nd);

    DET_TYPE hf_det = init_hf; 

    // using 1st state in CISD 
    calc_CISD(hf_det, b0, 2);
    b0.normalize();

    cout.precision(8);  
    
    sparse_solve_SI(b0,target_energy);

    //cout << "expectation val from dense vec:" << dense_b.transpose() * SpH * dense_b <<endl;; 
    return;

}
