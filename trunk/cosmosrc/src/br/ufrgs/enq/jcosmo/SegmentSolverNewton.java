package br.ufrgs.enq.jcosmo;

import no.uib.cipr.matrix.DenseMatrix;
import no.uib.cipr.matrix.DenseVector;
import no.uib.cipr.matrix.Vector;
import no.uib.cipr.matrix.Vector.Norm;

/**
 * Solve the SEGGAMMA by a full Newton method.
 * 
 * @param SEGGAMMA
 */
public class SegmentSolverNewton implements ISegmentSolver {
	
	int compseg;
	DenseMatrix jacobian;
	Vector residuals;
	
	public void solve(double PROFILE[], double factor, double SEGGAMMA[], double expDeltaW_RT[][], double tol){
		if(jacobian==null){
			compseg = PROFILE.length;
			jacobian = new DenseMatrix(compseg, compseg);
			residuals = new DenseVector(compseg);
		}
		
		// ITERATION FOR SEGMENT ACTIVITY COEF. (MIXTURE)
		
		// Solves the function vector
		// F_m = SEGGAMMA[m] * sum_n(PROFILE[n] * factor* SEGGAMMA[n] * expDeltaW_RT[m][n])
		int niter = 0;
		int maxiter = 20;

		while(true){
			for (int m = 0; m < compseg; m++) {
				double SUMMATION = 0.0;
				for(int n = 0; n < compseg; n++) {
					SUMMATION += PROFILE[n]*factor* SEGGAMMA[n] * expDeltaW_RT[m][n];
					
					// SEGGAMMA[m] * DSUMMATION/DSEGGAMMA[n]
					jacobian.set(m, n, SEGGAMMA[m] * (PROFILE[n]*factor*expDeltaW_RT[m][n]) );
				}
				jacobian.set(m, m, SUMMATION + jacobian.get(m, m));
				
				residuals.set(m, SEGGAMMA[m] * SUMMATION - 1.0);
			}
			
			double norm = residuals.norm(Norm.Two);
			if(norm <= tol || niter>maxiter)
				break;

			// update the vector accordingly to the Newton formula
			jacobian.solve(residuals, residuals);
			for (int i = 0; i < SEGGAMMA.length; i++)
				SEGGAMMA[i] -= residuals.get(i);

			++niter;
		}
		System.out.println("SEGGAMMA niter:" + niter);
	}
}
