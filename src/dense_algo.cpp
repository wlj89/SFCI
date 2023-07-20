/*
	Dense vector algorithms based on Eigen
	Comes in handy when testing new ideas
*/

#include "all.h"

using colvec = dynVec<DET_TYPE,FLOAT_TYPE>; 

void OffDiagGen_bit_e(const DET_TYPE, colvec&);
double DiagCal_bit(const DET_TYPE);


void comb_det(int n, int r, int *arr, int sz, vector<DET_TYPE>& basis)
{
    /*
		enumerate all determinants in <Sz>=0 block 
		
        sz = electronNum
    */
    DET_TYPE one = 1;

    for (int i = n; i >= r; i --)
    {
        arr[r - 1] = i;
        if (r > 1)
            comb_det(i - 1, r - 1, arr, sz, basis);
        else
        {
            DET_TYPE chi = 0; 

            int num = 0; 

            // <Sz>=0 
            for (int i = 0; i < sz; i ++)
            	if(arr[i]%2 == 0)
            		num ++;
            	else
            		num--;
                //chi[i] =arr[i]-1;
            
            if (num ==0)
            {
            	for (int i = 0; i < sz; i ++)
            	{
            		chi += ( one << (arr[i]-1) );    
            	}

                basis.push_back(chi);
			}

        }
    }
}

void generateBasis( vector<DET_TYPE>& basis, 
                    unordered_map<DET_TYPE,int>& basis_idx)
{
	/*
		generate basis 
	*/
    int n = orbitalNum; 
    int r = electronNum;
    int sz = r;
    int arr[electronNum];

    comb_det(n,r, arr, sz, basis);

    /*
        basis index table (determinant->index)
    */
    for (int i=0; i< basis.size(); i++ )
        basis_idx[basis[i]] = i; 
}


void fci(   vector<DET_TYPE>& basis,
            unordered_map<DET_TYPE,int>& basis_idx)
{
    /*
        FCI for small system
    */
    //void OffDiagGen_bit_e(const DET_TYPE, colvec&);
    //double DiagCal_bit(const DET_TYPE);

    int dim = basis.size();

    M HH(dim,dim);
    HH.setZero();

    M R(dim,dim);
    R.setIdentity();

    for(int i=0; i< dim; ++i)
    {
        colvec tgt;
        
        DET_TYPE current_det = basis[i]; 
        OffDiagGen_bit_e(current_det, tgt);
        tgt.merge(); 

        // diagonal
        HH(i,i) = DiagCal_bit(current_det);

        //off diag
        for (int j=0; j <tgt.len; j++)
        {
            int pos = basis_idx[tgt.basis[j]];
            HH(i,pos) = tgt.amp[j];
        }
    }

    cout.precision(12);

    SAES saes;
    saes.compute(HH);
    auto eigenval = saes.eigenvalues();
    auto eigenvec = saes.eigenvectors();

    cout <<"Ground state energy from FCI:"<<eigenval(0)<<endl;

    return;

}

void loadFromT( SpMat& SpH, 
                vector<DET_TYPE>& basis, 
                unordered_map<DET_TYPE,int>& basis_idx,
                double target_energy)
{
    /*
        fill the Row major sparse matrix 

        SpMat is set to be RowMajor! 
    */

    cout << "filling sparse matrix with target energy " << target_energy <<endl;

    unsigned int merge_bucket(dynVec<DET_TYPE,FLOAT_TYPE>&, unsigned int); 
    void introsort(dynVec<DET_TYPE,FLOAT_TYPE>&,int);


    int spaceDim = basis.size(); 
    int max_len = 0;  

    //#pragma omp parallel for
    // 21800 
    for (int i = 0; i< spaceDim; i++ )
    {
        if (i%10000==0) 
            cout << i << " rows loaded..."<<endl; 

        colvec tgt; 

        DET_TYPE current_det = basis[i]; 
        
        double current_diag_val = DiagCal_bit(current_det);
        SpH.insert(i,i) = (current_diag_val - target_energy);  

        OffDiagGen_bit_e(current_det,tgt);

        //cout << "len before merging:" << tgt.len <<endl; 
        //cout << "capacity of tgt:" << tgt.cap <<endl; 

        introsort(tgt,tgt.len);

        //tgt.checkOrder();

        unsigned int tmp_len = merge_bucket(tgt,tgt.len);
        
        tgt.len = tmp_len;
        
        //cout << tgt.len <<endl; 
        //cout << "offdiag done\n";
        
        //tgt.merge();

        for(int j=0; j< tgt.len; j++)
        {
            
            int pos = basis_idx[tgt.basis[j]];
            SpH.insert(i,pos) = tgt.amp[j];
        }
    }

    SpH.makeCompressed();   
    
    cout << "sparse matrix filled and compressed." <<endl;

}

void loadFromT_colmaj( SparseMatrix<double>& SpH, 
                vector<DET_TYPE>& basis, 
                unordered_map<DET_TYPE,int>& basis_idx)
{
    /*
        fill the col major sparse matrix 

        SpMat is set to be RowMajor! 
    */

    int spaceDim = basis.size(); 
    int max_len = 0; 

    //#pragma omp parallel for 
    for (int i =0; i< spaceDim; i++ )
    {
        if (i%10000==0) 
            cout << i << " rows filled..."<<endl; 

        colvec tgt; 

        DET_TYPE current_det = basis[i]; 
        
        double current_diag_val = DiagCal_bit(current_det);
        SpH.insert(i,i) = current_diag_val;  

        //cout << "dia done\n";

        OffDiagGen_bit_e(current_det,tgt);
        tgt.merge();

        for(int j=0; j< tgt.len; j++)
        {
            int pos = basis_idx[tgt.basis[j]];
            SpH.insert(i,pos) = tgt.amp[j];
        }
        
    }

    SpH.makeCompressed();   
    
    cout << "ColMajor sparse matrix filled and compressed." <<endl;

}

void krylov_solve_eigenval(SpMat& H, 
                          ColVec& v,
                              int spaceDim,
                              int num_iter,
                              int num_restart,
                              int num_states)
{
    /**
        Core of Lanzcos.
        Routine A(2,7) by C. C. Paige.
    */
    //cout << "num of threads used:" << nbThreads() <<endl<<endl;

    int numIterations = num_iter;
    int numRestarts = num_restart;

    int num_state_reached = 0; 
    double gamma;

    double energy_now; 
    double dE,energy_i,sigma_i;
    // stopping criteria 
    const double energy_threshold = 0.0000001;

    M KrylovBasis(spaceDim,numIterations);
    M KrylovM(numIterations,numIterations);

    KrylovM.setZero();
    KrylovBasis.setZero();
    
    M ritz_vec(spaceDim,numIterations);
    ritz_vec.setZero();
    /// set 0
    SAES saes;

    // record (i+1)th Ritz vector
    ColVec v_temp; 

    vector<double> lambda_list;
    vector<double> sigma_list;


    while (num_state_reached < num_states)
    {   

        /*
            Given the number of exterme eigenpairs to be sought. 


        */  
        energy_now = 0; 

        // num_restart is the max restart number 
        for(int res = 0; res < numRestarts; res++)
        {

            KrylovBasis.col(num_state_reached) = v;

            ColVec u = H* v;
            //cout <<"len u:"<<getNZ(u) <<endl; 

            ColVec w;
            ColVec v_next;
            
            double alpha = v.transpose() * u;
            KrylovM(num_state_reached,num_state_reached) = alpha;

            double beta;

            for (int i=1+num_state_reached; i<numIterations; i++)
            {
                /*
                    should re-orthogonalize
                */

                w = u - alpha * v;
                //cout <<"len w:"<<getNZ(w) <<endl; 
                cout << "gamma=" << gamma << endl;
                gamma = w.norm();
                //cout <<"gamma:" <<gamma<<endl;
                v_next = w / gamma;

                beta = gamma;

                u = H * v_next - beta * v;
                alpha = v_next.transpose() * u ;

                v = v_next; 
                
                /*
                for (int j = num_state_reached-1; j>=0; j--)
                    v = v - v.dot(ritz_vec.col(j)) * ritz_vec.col(j);

                for (int j = i-1; j>=0; j--)
                    v = v - v.dot(KrylovBasis.col(j))* KrylovBasis.col(j);

                v.normalize(); 
                */

                KrylovBasis.col(i) = v;

                KrylovM(i,i-1) = beta;
                KrylovM(i-1,i) = beta;
                KrylovM(i,i) = alpha;
            }

            saes.compute(KrylovM);
            auto lambda = saes.eigenvalues();
            auto eigenvec = saes.eigenvectors();
            

            //cout << "overlap matrix" <<endl;
            //cout << KrylovBasis.transpose() * KrylovBasis <<endl<<endl;;
            //cout << KrylovM <<endl;

            //cout << "basis overlap with preivous Ritz vec:" <<endl;
            //cout << ritz_vec.col(0).transpose() * KrylovBasis << endl;

            // orthginalize againt all previous Ritz vector
            
            v = KrylovBasis*eigenvec.col(num_state_reached);
            if (num_state_reached>0)
                for (int j = num_state_reached-1; j>=0; j--)
                    v = v - v.dot(ritz_vec.col(j)) * ritz_vec.col(j);
            v.normalize();

            //(i+1)th Ritz vector 
            v_temp = KrylovBasis*eigenvec.col(num_state_reached+1);
            if (num_state_reached>0)
                for(int j = num_state_reached-1; j>=0; j--)
                v_temp = v_temp - v_temp.dot(ritz_vec.col(j)) * ritz_vec.col(j);
            v_temp.normalize();

            //cout << "i & i+1 Ritz:" << v_temp.transpose()*v <<endl;;

            energy_i = lambda(num_state_reached);

            ColVec q = H * v - energy_i * v;
            sigma_i = q.norm();

            
            cout << "sigma:" << sigma_i<<endl;
            cout << "energy:" << energy_i<<endl;
            cout << "# round " << res<<endl;
            cout << "***************************************\n\n";

            dE = energy_now - energy_i; 
            energy_now = energy_i;

            if ( abs(dE) < energy_threshold)
            {
                cout << "energy discrepancy threshold reached. Terminating. " <<endl;
                break;
            }
        } 

        cout << "***************ATTENTION***************\n";
        cout << "calculation of state " << num_state_reached << " complete.\n";
        cout << "***************************************\n";

        lambda_list.push_back(energy_i); 
        sigma_list.push_back(sigma_i);

        KrylovM.setZero();
        KrylovBasis.setZero();

        // save the ith Ritz value
        for (int i=0; i< lambda_list.size(); i++)
            KrylovM(i,i) = lambda_list[i];
        
        // save the ith Ritz vector
        //KrylovBasis.col(num_state_reached) = v;
        ritz_vec.col(num_state_reached) = v; 
        //update v 
        /*
            enforcing orthogonalization of v0
        */
        v = v_temp; 
        
        v.normalize();

        //cout << "after enforcing orthogonalization" <<endl;
        //cout << KrylovBasis.col(num_state_reached).transpose() * v <<endl; 
        num_state_reached++;
    }

    cout << "energy of states\n";

    for(int i=0; i< num_states; i++)
    {
        cout << i << ":" << lambda_list[i]<< ", sigma:"<< sigma_list[i]<<endl; 
    }
}

void krylov_solve_linear(ColVec& x0,
                  ColVec& b,
                   SpMat& A,
                      int num_restart,
                      int num_iter/* dim of H matrix*/,
                   double residue_threshold)
{
    /*
        Krylov iterative solver of Ax = b 
        Restarted error projection 
    */

    M H(num_iter+1,num_iter); 
    M V(x0.size(),num_iter);

    H.setZero();
    V.setZero(); 

    double residue; 
    double r0_norm;

    for (int i = 0; i < num_restart; i++)
    {   

        cout << i << endl;

        ColVec r0 = b - A*x0;   
        ColVec v = r0.normalized(); 
        r0_norm = r0.norm(); 

        cout << "norm of r0: " << r0_norm <<endl; 

        V.col(0) = v;

        ColVec u = A*v;

        ColVec w;
        ColVec v_next;

        double alpha = v.transpose() * u; 

        H(0,0) = alpha;

        double beta, gamma;

        for (int j=1; j<num_iter+1; j++)
        {
            w = u - alpha* v;

            gamma = w.norm();

            v_next = w/gamma;

            beta = gamma;

            u = A*v_next - beta*v;

            alpha = v_next.transpose() * u;

            v= v_next;

            H(j,j-1) = beta;

            if (j < num_iter)
            {
                V.col(j) = v;
                
                H(j-1,j) = beta; 

                H(j,j) = alpha;
            }
        }   
        //cout << V.transpose() * V <<endl; 
        /*
            
            solve for ym = Hm^(-1)(beta,0,0,0...)^(T) -> error projection 
            
            this could be awkward since residue does not decrease monotonously

            should usw GMRES style least square 
        */

        //cout << H <<endl; 

        ColVec e1 = ColVec::Zero(num_iter+1);
        e1(0) = r0_norm;

        ColVec y = H.colPivHouseholderQr().solve(e1); 

        x0 = x0 + V*y; 

        //x0 = x0 + V * (H.inverse().col(0)*r0_norm);

        residue =  (b - A*x0).norm(); 

        cout << "after projection " << i+1 << ": " << residue <<endl; 

        if(residue < residue_threshold)
        {
            cout << "residue norm reached threshold. Terminating linear solver ." <<endl;
            break; 
        }
    }
    
    //cout << x0.norm() <<endl;
}



/*
void test_builtin_solver(   int spaceDim,
                            vector<DET_TYPE>& basis, 
                            unordered_map<DET_TYPE,int>& basis_idx)
{

    SparseMatrix<double> SpH(spaceDim,spaceDim);
    SpH.reserve(VectorXi::Constant(spaceDim,500));
    loadFromT_colmaj(SpH,basis,basis_idx);

    ColVec x(spaceDim);
    ColVec b = ColVec::Constant(spaceDim,0.1);

    SparseLU<SparseMatrix<double>, COLAMDOrdering<int>>  solver;

    solver.analyzePattern(SpH);

    solver.factorize(SpH); 

    x = solver.solve(b);

    cout << "residue:" <<(SpH*x - b).norm() <<endl;  
}

*/

void calc_diagonla_guess(ColVec& x, ColVec& b, SpMat& SpH)
{
    for (int i=0; i< b.size(); i++)
    {
        x(i) = b(i) / SpH.coeff(i,i);
    }
}


void solve_shift_and_invert(vector<DET_TYPE> basis,
                 unordered_map<DET_TYPE,int> basis_idx,
                                         int num_SI_iter,
                                         int num_solver_iter,
                                         int num_solver_restart,
                                      double target_energy,
                                      double residue_threshold)
{
    /*  
        shifted matrix     
    */
    int spaceDim = basis.size(); 
    int rowSpace = 2000;     
    SpMat SpH(spaceDim,spaceDim);
    SpMat SpH_original(spaceDim,spaceDim);

    SpH_original.reserve(VectorXi::Constant(spaceDim,rowSpace));
    SpH.reserve(VectorXi::Constant(spaceDim,rowSpace)); 


    loadFromT(SpH,basis,basis_idx,target_energy);

    SpH_original = SpH;

    for (int i=0; i <spaceDim; i++)
        SpH_original.coeffRef(i,i) += target_energy; 
    // initial guess of x. Should be more careful
    ColVec result;

    // initial guess of psi i.e. the final excited state 
    ColVec psi = ColVec::Random(spaceDim);     
    psi.normalize();

    for(int i = 0; i< num_SI_iter; i++)
    {

        //if (i==0)
        result = ColVec::Zero(spaceDim);

        //calc_diagonla_guess(result, psi, SpH);

        krylov_solve_linear(result,
                            psi,
                            SpH,
                            num_solver_restart,
                            num_solver_iter,
                            residue_threshold);   

        psi = result; 
        psi.normalize();

        ColVec HPsi = SpH_original * psi; 
        double energy_esti = (psi.transpose() * HPsi);
        //energy_esti /= psi.norm(); 

        //double sigma = (HPsi - energy_esti * psi).norm();    

        cout << "*****************************************\n";
        cout << "estimated energy:" << energy_esti <<endl;
        //cout << "sigma:" << sigma <<endl;   
        cout << "*****************************************\n";

    }

}

void dense_trunc(ColVec& x)
{
    int num_trunc = 0; 
    for(int i=0; i< x.size(); i++)
    {
        if (abs(x(i)) < 5e-5) 
        {
            x(i) = 0; 
            num_trunc ++; 
        }
    }

    cout << "size of the truncated: " << x.size() - num_trunc <<endl; 
}

void krylov_solve_linear_nonorthog( ColVec& x0,
                                    ColVec& b,
                                    SpMat& SpH)
{
    /*
        unorthogonal basis Krylov linear solver
    */    
    int spaceDim = x0.size(); 

    x0 = ColVec::Zero(x0.size());
    ColVec w(x0.size());
    //x0(spaceDim/22+123) = 1.0;

    //ColVec w(x0.size());

    x0.normalize();

    ColVec r0 = b - SpH * x0;

    double r0_norm = r0.norm(); 

    cout << "initial residue:" << r0_norm <<endl; 
    int num_iter = 15; 
    M hessenberg(num_iter+1, num_iter);
    M krylov_basis(spaceDim,num_iter);
    M S(num_iter+1,num_iter+1);

    hessenberg.setZero();

    krylov_basis.setZero();
    krylov_basis.col(0) = r0;

    SAES saes;

    double w_norm_final;

    S(0,0) = r0.norm();

    cout << "beep\n" <<endl; 

    for(int i=0; i<num_iter; i++)
    {
        w = SpH * krylov_basis.col(i); 

        vector<double> Hij_temp(num_iter,0);
        // alpha_i 
        //hessenberg(i,i) = w.dot(krylov_basis.col(i));

        for (int j=0; j<=i; j++)
        {
            Hij_temp[j] = w.dot(krylov_basis.col(j));
        }

        //for (int i)
        /*
            full orthogonalization
        */
        /*
        for (int j=0; j<=i; j++)
            w = w - Hij_temp[j] * krylov_basis.col(j);
        */
        /*
            tri-orthogonalization
        */
        
        if (i>0)    
        {
            w = w - Hij_temp[i] * krylov_basis.col(i)  - 
                    Hij_temp[i-1] * krylov_basis.col(i-1);
        }   
        else
        {
            w = w - Hij_temp[i] * krylov_basis.col(i);
        }   
        
        // apply truncation
        cout << i <<endl;  
        //dense_trunc(w);

        for (int j=i; j>=0; j--)
        {
            hessenberg(i,j) = Hij_temp[j];
            hessenberg(j,i) = hessenberg(i,j); 
        }

        dense_trunc(w);
        w.normalize(); 
            

        if (i < (num_iter-1))
        {      
            krylov_basis.col(i+1) = w; 
        }
        

        // compute S 

        S(i+1,i+1) = w.norm();

        for (int j = 0; j<=i; j++)
        {
            S(i+1,j) = krylov_basis.col(j).dot(w);
            S(j,i+1) = S(i+1,j);
        }
    }  

    ColVec w_final = SpH* w;

    for (int i=0; i< num_iter; i++)
    {
        hessenberg(num_iter,i) = w_final.dot(krylov_basis.col(i));
    }

    //auto hess_tilda = hessenberg.block(0,0,num_iter,num_iter-1);
    /*
    ColVec e1 = ColVec::Zero(num_iter+1);
    
    e1(0) = r0_norm;  
    */
    
    ColVec e1 = S.col(0);
    e1 *= r0_norm;
    
    ColVec y =   (S*hessenberg).colPivHouseholderQr().solve(e1); 
    
    for (int i=0; i< num_iter; i++)
        x0 = x0 + krylov_basis.col(i)*y(i); 
    //x0 = x0 + V * (H.inverse().col(0)*r0_norm);

    double residue =  (b - SpH*x0).norm(); 

    cout << "truncating resultant x0" <<endl;
    dense_trunc(x0);

    double residue_trunc = (b - SpH*x0).norm(); 

    cout << "residue: "<<residue <<endl; 
    cout << "residue trunc: " << residue_trunc <<endl;


    cout << "overlap\n" << S <<endl;

    cout << "hessenberg\n" << hessenberg <<endl;
    /*
    cout << "calculated hessenberg:\n"; 
    cout << hessenberg <<endl; 
    cout << "actual hessenberg:\n";
    cout << krylov_basis.transpose() * SpH * krylov_basis <<endl; 
    */
    /*
    cout << "calculated residue:" << w_norm_final * eigenvec.col(0)[num_iter-1]<<endl;;
    cout << "calculated expectation:" << lambda(0)<<endl; 
    ColVec result = krylov_basis * eigenvec.col(0); 

    cout << "solution's norm:" << result.norm() <<endl; 
    ColVec temp = SpH * result; 

    cout << "actual expectation:" << result.dot(temp) <<endl;

    cout << "actual residue:" << sqrt(temp.dot(temp) - lambda(0)*lambda(0))<<endl;
    //cout << lambda(0) <<endl; 
    //cout << "Hessenberg matrox:" << endl;
    //cout << hessenberg <<endl;  

    cout << "overlap\n";
    cout << krylov_basis.transpose() * krylov_basis << endl;
    */

}


void test(SpMat& SpH, int spaceDim)
{   
    /*
        a more efficient projection routine
    */

    //looks like a HF solution
    ColVec x0 = ColVec::Zero(spaceDim);
    x0(spaceDim/22+123) = 1.0;

    ColVec w(spaceDim);

    x0.normalize();

    int num_iter = 10; 
    M hessenberg(num_iter, num_iter);
    M krylov_basis(spaceDim,num_iter);
    
    hessenberg.setZero();

    krylov_basis.setZero();
    krylov_basis.col(0) = x0;

    SAES saes;

    double w_norm_final;

    for(int i=0; i<num_iter; i++)
    {
        w = SpH * krylov_basis.col(i); 

        // alpha_i 
        hessenberg(i,i) = w.dot(krylov_basis.col(i));
        
        if (i>0)    
        {
            // A*v_i - alpha_i *v_i - beta_{i-1} * v_{i-1}
            w = w - hessenberg(i,i) * krylov_basis.col(i)  - 
                    hessenberg(i,i-1) * krylov_basis.col(i-1);
        }   
        else
        {
            w = w - hessenberg(i,i) * krylov_basis.col(i);
        }

        // test non-orthogonal basis

        if (i < (num_iter-1)) 
        {   
            hessenberg(i,i+1) = w.norm(); 
            hessenberg(i+1,i) = hessenberg(i,i+1);

            krylov_basis.col(i+1) = w/hessenberg(i,i+1);
        }
        else
        {
            w_norm_final = w.norm(); 
        }
    }

    saes.compute(hessenberg);

    auto lambda = saes.eigenvalues();
    auto eigenvec = saes.eigenvectors();

    /*
    cout << "calculated hessenberg:\n"; 
    cout << hessenberg <<endl; 
    cout << "actual hessenberg:\n";
    cout << krylov_basis.transpose() * SpH * krylov_basis <<endl; 
    */

    cout << "calculated residue:" << w_norm_final * eigenvec.col(0)[num_iter-1]<<endl;;
    cout << "calculated expectation:" << lambda(0)<<endl; 
    ColVec result = krylov_basis * eigenvec.col(0); 

    cout << "solution's norm:" << result.norm() <<endl; 
    ColVec temp = SpH * result; 

    cout << "actual expectation:" << result.dot(temp) <<endl;

    cout << "actual residue:" << sqrt(temp.dot(temp) - lambda(0)*lambda(0))<<endl;
    //cout << lambda(0) <<endl; 
    //cout << "Hessenberg matrox:" << endl;
    //cout << hessenberg <<endl;  

    cout << "overlap\n";
    cout << krylov_basis.transpose() * krylov_basis << endl;
}

void dense_main()
{
    void sparse_solve_SI(int, vector<DET_TYPE>&,colvec&);

    setNbThreads(4);

    unsigned int seed = 123123; 
    //srand((unsigned int) time(0));
    srand(seed);

    cout.precision(7);  

    vector<DET_TYPE> basis;    
    unordered_map<DET_TYPE,int> basis_idx;

    generateBasis(basis,basis_idx);
    cout << "space size:" << " "<< basis.size() << endl;

    int spaceDim = basis.size();

    //fci(basis, basis_idx); 

    //test_builtin_solver(spaceDim,basis,basis_idx);

    double res_threshold = 1E-6; 
    double target_energy = -74.71;
    int num_SI_iter = 20; 
    int num_solver_iter = 20;
    int num_solver_restart = 10; 
    
    /*
    solve_shift_and_invert( basis,
                            basis_idx,
                            num_SI_iter,
                            num_solver_iter,
                            num_solver_restart,
                            target_energy,
                            res_threshold);
    return; 
    */

    int rowSpace = 2000;     
    SpMat SpH(spaceDim,spaceDim);
    SpH.reserve(VectorXi::Constant(spaceDim,rowSpace)); 
    loadFromT(SpH,basis,basis_idx,target_energy);

    //ColVec b = ColVec::Random(spaceDim);
    
    ColVec x0(spaceDim);
    /*
    x0 = ColVec::Random(spaceDim);
    x0.normalize();
    b.normalize(); 
    */
    ColVec b = ColVec::Zero(spaceDim);
    b(spaceDim/223+61) = 1.0;

    /*take the diagonal guess*/
    //calc_diagonla_guess(x0,b,SpH);
    krylov_solve_linear_nonorthog(x0,b,SpH);
    //test(SpH,spaceDim); 

    return;
    /*
    krylov_solve_eigenval( SpH,
                           x0,
                           spaceDim,
                           20, //dim of the tridiagonal matrix
                           80,
                           6);
    */
    //parallel_SpMVM_testse(SpH, x0);

    //krylov_solve_linear(x0,b,SpH,10,100,res_threshold);
    
}

