package br.ufrgs.enq.jcosmo;

import org.netlib.blas.BLAS;

/**
 * Solve the SEGGAMMA by successive substitutions.
 * 
 * @param SEGGAMMA
 */
public class SegmentSolverSimple implements ISegmentSolver {
	
	public void solve(double PROFILE[], double factor, double SEGGAMMA[], double expDeltaW_RT[][], double tol){
		int nsegments = PROFILE.length;
		
		int niter = 0;
		int maxiter = 100;
		double norm = -1;
		BLAS blas = BLAS.getInstance();
		
		while(true){
			for (int m = 0; m < nsegments; m++) {
//			for (int m = compseg-1; m >= 0; m--) {
				double SUMMATION = 0.0;
				for(int n = 0; n < nsegments; n++) {
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
//		System.out.println("SEGGAMMA niter:" + niter);
	}
	
	
	public void solveMulti(double PROFILE[][], double factor, double SEGGAMMA[][], double expDeltaW_RT[][][], double tol){
		int ndescriptors = PROFILE.length;
		int nsegments = PROFILE[0].length;
		
		int niter = 0;
		int maxiter = 100;
		double norm = -1;
		double normVector[] = new double[ndescriptors];
		BLAS blas = BLAS.getInstance();
		
		while(true){
			for (int d = 0; d < ndescriptors; d++) {
				for (int m = 0; m < nsegments; m++) {
					double SUMMATION = 0.0;
					for(int n = 0; n < nsegments; n++) {
						SUMMATION += PROFILE[d][n]*factor* SEGGAMMA[d][n] * expDeltaW_RT[d][m][n];
					}
					// substitute with the new value
					SEGGAMMA[d][m] = 1.0/SUMMATION;
				}
				normVector[d] = blas.dnrm2(nsegments, SEGGAMMA[d], 1);
			}
			
			double newnorm = blas.dnrm2(ndescriptors, normVector, 1);
			if(Math.abs((norm - newnorm)/newnorm) <= tol || niter>maxiter)
				break;
			norm = newnorm;

			++niter;
		}
//		System.out.println("SEGGAMMA niter:" + niter);
	}

}
