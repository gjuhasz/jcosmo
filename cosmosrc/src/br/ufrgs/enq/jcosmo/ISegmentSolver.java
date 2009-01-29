package br.ufrgs.enq.jcosmo;

/**
 * Interface to be implemented by segment activity solvers.
 * 
 * <p>The function to be solved has the following form:
 * <p> F_m = SEGGAMMA[m] * sum_n(PROFILE[n] * factor* SEGGAMMA[n] * expDeltaW_RT[m][n])
 * 
 */
public interface ISegmentSolver {
	
	/**
	 * Solve the segment activity function.
	 * 
	 * @param PROFILE
	 * @param factor
	 * @param SEGGAMMA the segment activity vector (should contain the solution on exit)
	 * @param expDeltaW_RT
	 * @param tol
	 */
	public void solve(double PROFILE[], double factor, double SEGGAMMA[], double expDeltaW_RT[][], double tol);

}