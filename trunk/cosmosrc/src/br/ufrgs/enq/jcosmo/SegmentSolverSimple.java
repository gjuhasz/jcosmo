package br.ufrgs.enq.jcosmo;

import org.netlib.blas.BLAS;

/**
 * Solve the SEGGAMMA by successive substitutions.
 * 
 * @param SEGGAMMA
 */
public class SegmentSolverSimple implements ISegmentSolver {
	
	int compseg;

	public void solve(double PROFILE[], double factor, double SEGGAMMA[], double expDeltaW_RT[][], double tol){
		compseg = PROFILE.length;
		
		int niter = 0;
		int maxiter = 100;
		double norm = -1;
		BLAS blas = BLAS.getInstance();
		
		while(true){
			for (int m = 0; m < compseg; m++) {
				double SUMMATION = 0.0;
				for(int n = 0; n < compseg; n++) {
					SUMMATION += PROFILE[n]*factor* SEGGAMMA[n] * expDeltaW_RT[m][n];
				}
				// substitute with the new value
				SEGGAMMA[m] = 1.0/SUMMATION;
			}
			
			double newnorm = blas.dnrm2(SEGGAMMA.length, SEGGAMMA, 1);
			if(Math.abs((norm - newnorm)/newnorm) <= tol || niter>maxiter)
				break;
			norm = newnorm;

			++niter;
		}
		System.out.println("SEGGAMMA niter:" + niter);
	}
}
