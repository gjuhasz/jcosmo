package br.eng.rps.cosmo;


/**
 * COSMO-SAC activity model.
 * 
 * This activity implementation uses the sigma profiles of pure components to calculate
 * the liquid-phase activity coefficients in a solution.
 * <p>
 * This code is based on code available at http://www.design.che.vt.edu/COSMO/
 * <p>
 * There are several improvements on the code:
 * <ul>
 * <li>Code was translated for Java
 * <li>It can evaluate infinite dilution (x=0) activity coeficients (original code
 * had a bug and returned NaN for this case)
 * <li>There are some changes on the segment activity converging mechanism (much better performance)
 * </ul>
 * 
 * @author Rafael de Pelegrini Soares
 * 
 */
public class COSMOSAC {

	static final double EO = 2.395e-4;
	public static final double AEFFPRIME = 7.5;
	double aEffPrime = AEFFPRIME;
	
	static final double RGAS = 0.001987;
	public static final double VNORM = 66.69;
	public static final double ANORM = 79.53;

	double vnorm = VNORM;
	double anorm = ANORM;

	// (LIN AND SANDLER USE A CONSTANT FPOL WHICH YEILDS EPS=3.68)
	public static final double EPSILON = 3.667;
	
	/** Default coordination number (KLAMT USED 7.2) */
	public static final double COORD = 10.0;
	double coord = COORD;
	
	public static final double SIGMAHB = 0.0084;
	public static final double CHB = 85580.0;
	
	double sigmaHB = SIGMAHB;
	double cHB = CHB;

	double alphaPrime;
	
	private static final double TOLERANCE = 1e-6;

	protected double T;
	protected double inv_RT;
	protected double[] z;

	// The number of components
	int ncomps;
	int compseg;

	/// Cavity volume for each substance
	double [] VCOSMO;
	double [] ACOSMO;

	double [] RNORM;
	double [] QNORM;
	double [] li;

	double [] PROFILE;

	double [][] sigma;
	double [] charge;
	double deltaW[][];
	double SEGGAMMA[];
	double SEGGAMMAPR[][];
	
	/** Flag if the Staverman-Guggenheim term is ignored */
	boolean ignoreSG = false;

	/**
	 * @see #setParameters(double[], double[], double[][])
	 */
	public COSMOSAC(){
	}
	
	
	public double getAEffPrime() {
		return aEffPrime;
	}


	public void setAEffPrime(double effPrime) {
		aEffPrime = effPrime;
	}


	public double getCoord() {
		return coord;
	}


	public void setCoord(double coord) {
		this.coord = coord;
	}

	/**
	 * Creates a new COSMO-SAC object given the cavity volumes and sigma profiles.
	 * 
	 * Multi-compound mixtures are supported.
	 * 
	 * If {@link #charge} has length nseg then {@link #sigma} should have dimension
	 * of [ncomps][nseg], were ncomps is the number of components in mixture.
	 * 
	 * @param cavityVolume the cavity volume for each compound, length [ncomp]
	 * @param charge the charge density, length [nseg]
	 * @param sigma the scaled sigma profile, dimensions [ncomp][nseg]
	 */
	public void setParameters(double []cavityVolume, double []charge, double [][] sigma){
		assert cavityVolume.length == sigma.length;

		this.ncomps = cavityVolume.length;
		this.compseg = charge.length;
		this.VCOSMO = cavityVolume;
		this.charge = charge;
		this.sigma = sigma;

		this.T = 300;

		z = new double[ncomps];
		for (int i = 0; i < ncomps; i++)
			z[i] = 1.0/ncomps;

		rebuildInternals();
	}

	/**
	 * Sets the system temperature in Kelvin
	 * @param T the temperature
	 */
	public void setTemperature(double T){
		this.T = T;
		this.inv_RT = 1.0 / (RGAS * T);
	}
	
	/**
	 * Sets the composition of the mixture in molar basis.
	 * @param z the molar composition.
	 */
	public void setComposition(double []z){
		for (int i = 0; i < ncomps; i++)
			this.z[i] = z[i];
		
		// CALCULATE THE MIXTURE SIGMA PROFILE
		double denom = 0;
		for (int i = 0; i < ncomps; i++) {
			denom += z[i]*ACOSMO[i];
		}
		for(int m=0; m<compseg; ++m){
			double numer = 0;
			for (int i = 0; i < ncomps; i++) {
				numer += z[i]*sigma[i][m]; 
			}
			PROFILE[m] = numer/denom;
		}
	}

	/**
	 * Reallocates all mixture dependent variables.
	 * This function calculates all quantities which do not depend on composition
	 * or temperature.
	 */
	private void rebuildInternals(){
		double alpha = 0.3*Math.pow(aEffPrime, 1.5)/EO;
		double fpol = (EPSILON-1.0)/(EPSILON+0.5);
		alphaPrime = fpol*alpha;
		
		ACOSMO = new double[ncomps];
		PROFILE = new double[compseg];
		deltaW = new double[compseg][];
		for(int i=0; i<compseg; ++i){
			deltaW[i] = new double[compseg];
		}

		SEGGAMMA = new double[compseg];
		SEGGAMMAPR = new double[ncomps][];
		RNORM = new double[ncomps];
		QNORM = new double[ncomps];
		li = new double[ncomps];
		for(int i=0; i<ncomps; ++i){
			SEGGAMMAPR[i] = new double[compseg];
			ACOSMO[i] = 0;
			for (int m = 0; m < sigma[i].length; m++) {
				ACOSMO[i] += sigma[i][m];
				
				// initialize all SEGGAMMAPR (reused between calculations)
				SEGGAMMAPR[i][m] = 1.0;
			}
		}
		
		double SIGMAACC = 0, SIGMADON = 0;
		double chargemn = 0;
		for(int m=0; m<compseg; ++m){
			// initialize all SEGGAMMA (reused between calculations)
			SEGGAMMA[m] = 1.0;

			for(int n=0; n<compseg; ++n){
				if(charge[m]>=charge[n]){
					SIGMAACC = charge[m];
					SIGMADON = charge[n];
				}
				else {
					SIGMADON = charge[m];
					SIGMAACC = charge[n];
				}
				chargemn = charge[m]+charge[n];
				deltaW[m][n] = (alphaPrime/2.0)*chargemn*chargemn +
					cHB * Math.max(0.0, SIGMAACC - sigmaHB)*Math.min(0.0, SIGMADON + sigmaHB);
			}
		}
		
		// some more composition independent properties
		for(int i=0; i<ncomps; ++i){
			RNORM[i] = VCOSMO[i]/vnorm;
			QNORM[i] = ACOSMO[i]/anorm;
			li[i] = (coord/2.0)*(RNORM[i]-QNORM[i])-(RNORM[i]-1.0);
		}
	}
	
	/**
	 * Calculates the activity coefficients.
	 * 
	 * @param lnGama the vector to put the results in.
	 */
	public void activityCoefficient(double[] lnGama){
		activityCoefficientLn(lnGama, 0);
	}

	/**
	 * Calculates the activity coefficients.
	 * 
	 * @param lnGama the vector to put the results in.
	 */
	public void activityCoefficientLn(double[] lnGama, int start){
		// ITERATION FOR SEGMENT ACTIVITY COEF. (MIXTURE)
		double norm = -1;
		int niter = 0;
		while(true){ 
			for (int m = 0; m < compseg; m++) {
				double SUMMATION = 0.0;
				for(int n = 0; n < compseg; n++) {
					SUMMATION += PROFILE[n]* SEGGAMMA[n] * Math.exp(-deltaW[m][n] * inv_RT);
				}
				// NOTE: starting value comes from last calculation
				// update the solution
				SEGGAMMA[m]=1.0/SUMMATION;
			}
			double newnorm = blas_dnrm2(SEGGAMMA.length, SEGGAMMA, 1);
			if(Math.abs((norm - newnorm)/newnorm) <= TOLERANCE)
				break;
			norm = newnorm;
			++niter;
		}
//		System.out.println("SEGGAMMA niter:" + niter + " norm:" + norm);


		// ITERATION FOR SEGMENT ACITIVITY COEF (PURE SPECIES)
		for(int i=0; i<ncomps; ++i){
			niter = 0;
			while(true){
				++niter;
				for(int m=0; m<compseg; ++m){
					double SUMMATION = 0.0;
					for(int n=0; n<compseg; ++n){
						SUMMATION += sigma[i][n]/ACOSMO[i]*SEGGAMMAPR[i][n]*Math.exp(-deltaW[m][n] * inv_RT);
					}
					// NOTE: starting value comes from last calculation
					// update with a new value
					SEGGAMMAPR[i][m] = 1.0/SUMMATION;
				}
				double newnorm = blas_dnrm2(SEGGAMMAPR[i].length, SEGGAMMAPR[i], 1);
				if(Math.abs((norm - newnorm)/newnorm) <= TOLERANCE)
					break;
				norm = newnorm;
			}
//			System.out.println("SEGGAMMA I:" + i + " niter:" + niter);
		}

		// THE STAVERMAN-GUGGENHEIM EQUATION
		double BOTTHETA = 0;
		double BOTPHI = 0;
		for(int i=0; i<ncomps; ++i){
			BOTTHETA += z[i]*QNORM[i];
			BOTPHI += z[i]*RNORM[i];
		}

		double sum_zi_li = 0;
		for(int i=0; i<ncomps; ++i)
			sum_zi_li += z[i]*li[i];

		for(int i=0; i<ncomps; ++i){
			double phi_z = RNORM[i]/BOTPHI;
			
			double theta_phi = QNORM[i]/BOTTHETA/(RNORM[i]/BOTPHI);

			// GAMMASGI IS ACTUALLY LNGAMMASG
			double lnGammaSG = 0;
			if(!ignoreSG){
				lnGammaSG += Math.log(phi_z)+
				(coord/2)*QNORM[i]*Math.log(theta_phi)
				+ li[i] - (phi_z)*sum_zi_li;
			}

			// CALCULATION OF GAMMAS
			double lnGammaRestoration = 0.0;
			for(int m=0; m<compseg; ++m){
				lnGammaRestoration += (sigma[i][m]/aEffPrime)*(Math.log(SEGGAMMA[m]/(SEGGAMMAPR[i][m])));
			}
			lnGama[i+start] = lnGammaRestoration + lnGammaSG;
		}
	}
	
	/**
	 * @param n the input vector length
	 * @param x the input vector
	 * @param incx the increment (usually 1)
	 * @return the euclidean norm of the given vector
	 */
	private static double blas_dnrm2(int n, double[] x, int incx){
		double sum = 0;
		for (int i = 0; i < n; i+= incx)
			sum += x[i] * x[i];
		return Math.sqrt(sum);
	}

	public boolean isIgnoreSG() {
		return ignoreSG;
	}

	/**
	 * Sets if the Staverman-Guggenheim term should be ignored.
	 * @param ignoreSG <code>true</code> if the term should be ignored
	 */
	public void setIgnoreSG(boolean ignoreSG) {
		this.ignoreSG = ignoreSG;
	}


	public double getAnorm() {
		return anorm;
	}


	public void setAnorm(double anorm) {
		this.anorm = anorm;
	}


	public double getVnorm() {
		return vnorm;
	}


	public void setVnorm(double vnorm) {
		this.vnorm = vnorm;
	}


	public double getCHB() {
		return cHB;
	}


	public void setCHB(double chb) {
		cHB = chb;
	}


	public double getSigmaHB() {
		return sigmaHB;
	}


	public void setSigmaHB(double sigmaHB) {
		this.sigmaHB = sigmaHB;
	}
}
