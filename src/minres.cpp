#include "all.h"

// define tolerance in head file

using colvec = dynVec<DET_TYPE, FLOAT_TYPE>;

void calExcitedEnergy(colvec& start_vec, 
                          int restartNum, 
                          int iterNum, 
                   FLOAT_TYPE tolerance,
                   FLOAT_TYPE lambda)
    
{
    int N = Nd; 

    colvec psi(N);
    colvec psi_guess(N);
    colvec psi_result(N);

    psi.naiveCopy(start_vec);
    psi_guess.setZeroVec(); 
    // randimzed psi 

    FLOAT_TYPE energy = calEnergy(psi); 

    cout << "initial energy: " << energy <<endl; 

    auto start = system_clock::now();
    auto end = system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end-start;

    for(int i=0; i< restartNum; i++)
    {
        
        cout << "#############################################\n";
		cout << "iteration " << i <<endl; 
        //    At some point, use the result as the guess for next iteration.
        start = system_clock::now(); 

        minres_sp(psi,psi_guess,psi_result,iterNum,tolerance,lambda); 
        psi_result.normalize(); 

        psi.naiveCopy(psi_result);
        
		//psi_guess.setZeroVec(); 
		psi_guess.naiveCopy(psi_result);

        energy = calEnergy(psi);

        end = system_clock::now();
        elapsed_seconds = end-start;
        cout << "time of one iteration: " <<elapsed_seconds.count()<<endl; 

        cout << "energy:  "<< energy <<endl;

        cout << "#############################################\n\n";
    }

    /// write wavefunction
    string output = "wf_self_as_guess.out";
    ofstream file; 
    
    file.open(output, ios::trunc);
    
    file << psi_result.len <<endl; 
    
    for (int i=0; i<psi_result.len; i++)
        file <<psi_result.basis[i]<<' '<< std::setprecision(12) << psi_result.amp[i] <<endl; 
    
    file.close(); 
}

void minres_sp(dynVec<DET_TYPE, FLOAT_TYPE>& rhs, 
               dynVec<DET_TYPE, FLOAT_TYPE>& guess, 
               dynVec<DET_TYPE, FLOAT_TYPE>& result, 
                                         int iters, 
                                  FLOAT_TYPE tolerance,
                                  FLOAT_TYPE lambda)
{
	/*
        Original code credit to Dr. Giamoco Po @ Miami U 
        MINRES with Spv and SpM  
        
        Note: in the following implementation, guess vector must be zero vector. 
        If one wants to reuse previous result, the expression of rhsNorm2 should be changed accordingly. 
        
        (For more details, refer to Chap 6.5.3 of Iterative Method for Sparse Linear System by Yousef Saad)

        For Dot product, normalization, addition/substraction .etc. 
        This is also the reserved mem for any SpV, though the actual usage coud be less. 
             Nd = 10000000 ? 5000000
        
        For SpVSpM:
            nDet_spvspm = 1000000  
        
        Each SpV must have strictly descending bit-string at ANY TIME.
        Bit string is more effcient than the amplitude ordered 
        
    */

    //parameters set-up 
    unsigned int maxIters = iters;  // initialize maxIters to iters
    unsigned int N = Nd;    // the size of the matrix

	FLOAT_TYPE diag;
    FLOAT_TYPE tol_error = tolerance;
    FLOAT_TYPE rhsNorm2 = rhs.norm(); 
    rhsNorm2 = rhsNorm2 * rhsNorm2; 
    FLOAT_TYPE threshold2 = tol_error*tol_error*rhsNorm2; // convergence threshold (compared to residualNorm2)
	
    //cout <<threshold2 <<endl; 	
    //solution

    colvec x(N);
    colvec v_old(N);    
    colvec v(N);
    colvec v_new(N); 
    colvec w(N);
    colvec w_new(N);

    x.naiveCopy(guess); 
    
    //ColVec v_old(N); // will be initialized inside loop
    //v = ColVec::Zero(N); //initialize 

    /*
    ColVec v_new = rhs-mat*x; //initialize v_new    
    */    
    SpMVMtply_bucket(x,v_new, lambda,-1.0,diag);
    //cout << "asd" <<endl;
	v_new.addTwo(rhs,1.0); 
	//cout << "asd1" <<endl; 
    FLOAT_TYPE residualNorm2 = pow(v_new.norm(),2);
    w_new.naiveCopy(v_new);
    
    FLOAT_TYPE beta_new2 = v_new.dot(w_new);
    FLOAT_TYPE beta_new = sqrt(beta_new2);
    FLOAT_TYPE beta_one = beta_new;

    v_new.scalarMtply(1.0/beta_new);
    w_new.scalarMtply(1.0/beta_new);

    FLOAT_TYPE c = 1.0; // the cosine of the Givens rotation
    FLOAT_TYPE c_old = 1.0;
    FLOAT_TYPE s = 0.0; // the sine of the Givens rotation
    FLOAT_TYPE s_old = 0.0; // the sine of the Givens rotation
    FLOAT_TYPE eta = 1.0;
	//FLOAT_TYPE diag; 

    colvec p_oold(N); // will be initialized in loop
    colvec p_old(N);
    colvec p(N);
    
    p_old.setZeroVec();
    p.setZeroVec(); 
	
	//return ; 
    
	int IterNum = 0; 

    while ( IterNum < maxIters )
    {
        FLOAT_TYPE beta = beta_new;
        
        v_old.naiveCopy(v);
        v.naiveCopy(v_new);
        w.naiveCopy(w_new);

        /*
        v_new = mat*w - beta*v_old; // compute v_new
        */
		
        SpMVMtply_bucket(w,v_new,lambda,1.0,diag);
		
	    //cout << beta<<endl;	
		v_new.addTwo(v_old,-beta);
		//cout  << "beep" <<endl; 
        
		FLOAT_TYPE alpha = v_new.dot(w);
         
        v_new.addTwo(v,-alpha);
        w_new.naiveCopy(v_new); 
		
		 
        beta_new2 = v_new.dot(w_new); // compute beta_new
        beta_new = sqrt(beta_new2);

        v_new.scalarMtply(1.0/beta_new);
        w_new.scalarMtply(1.0/beta_new);
		
		//cout << "beep1" <<endl;
        const double r2 =s*alpha+c*c_old*beta; // s, s_old, c and c_old are still from previous iteration
        const double r3 =s_old*beta; // s, s_old, c and c_old are still from previous iteration
        const double r1_hat=c*alpha-c_old*s*beta;
        const double r1 =sqrt( pow(r1_hat,2) + pow(beta_new,2) );
        c_old = c; // store for next iteration
        s_old = s; // store for next iteration
        c = r1_hat/r1; // new cosine  
        s = beta_new/r1; // new sine

        p_oold.naiveCopy(p_old);
        p_old.naiveCopy(p);
        p.addThree(w,p_old,p_oold,1.0/r1,-r2/r1,-r3/r1);
		//cout << "beep2"<<endl;
        x.addTwo(p,beta_one*c*eta);

        residualNorm2 *= s*s;   
		cout << "res norm:" << residualNorm2 <<endl;//<<endl; 
		cout << IterNum << " is done\n";
		cout << "*****************************************\n";
        if ( residualNorm2 < tol_error)
            break;  
        
        eta=-s*eta; // update eta
        IterNum++; // increment iteration number (for output purposes)
    }   

    result.naiveCopy(x); 

    //double actual_error = sqrt(residualNorm2/rhsNorm2);
    //cout << "error:" <<actual_error<<endl;
    //cout << "iter num:" << IterNum<<endl; 

    return;
}

