/*
 * Copyright 2008 Rafael de Pelegrini Soares and Renan Pereira Gerber
 * 
 * This file is part of JCosmo.
 * 
 * JCosmo is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * JCosmo is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 * 
 * You should have received a copy of the GNU Lesser General Public License
 * along with JCosmo.  If not, see <http://www.gnu.org/licenses/>.
 */

package br.ufrgs.enq.jcosmo;




/**
 * COSMO-SAC activity model with multiple descriptors.
 * 
 * This activity implementation uses multiple domain sigma profiles
 * of pure components to calculate the liquid-phase activity coefficients in a solution.
 * <p>
 * This code was initially based on code available at http://www.design.che.vt.edu/COSMO/
 * but currently only a very small portion of the original code remains.
 * 
 * @author Rafael de Pelegrini Soares
 * 
 */
public class COSMOSACMulti {

	static final double E0 = 2.395e-4;
	public static final double FPOL = 0.6917; // Ind. Eng. Chem. Res., Vol. 46, No. 22, 2007
	/** Default averaging radius */
	public static final double RAV = 0.81764200000000;

	double rav = RAV;
	double [][]fpol;
	double []beta;

	static final double RGAS = 0.001987; // kcal/mol/K
	public static final double VNORM = 66.69;
	public static final double ANORM = 79.53;

	double vnorm = VNORM;
	double anorm = ANORM;
	double rPower = 1;

	/** Default coordination number (KLAMT USED 7.2) */
	public static final double COORD = 10.0;
	double coord = COORD;

	public static final double SIGMAHB = 0.0084;
	public static final double CHB = 85580.0;
	
	public static final int HBTYPE_ORIGINAL = 0;
	public static final int HBTYPE_HSIEH = 1;
	public static final int HBTYPE_FIXED = 2;
	public static final int HBTYPE_COVALENT = 3;
	
	int hbType = HBTYPE_ORIGINAL;
	double sigmaHB = SIGMAHB;
	double sigmaHB2 = SIGMAHB;
	double sigmaHB3 = 1;
	double [][]cHB;
	
	private static double SIGMA_BOUND = 0.025;
	
	double alpha;

	private static final double TOLERANCE = 1e-4;

	protected double T;
	protected double inv_RT;
	protected double[] z;

	/** The number of compounds in the mixture */
	int ncomps;
	/** The number of independent segments to be considered by the model */ 
	int nsegments;
	/** The number of descriptors to be considered by the model */ 
	int ndescriptors;

	/// Cavity volume for each substance
	double [] VCOSMO;
	double [] ACOSMO;

	double [] RNORM;
	double [] QNORM;
	double [] li;

	double [][] PROFILE;

	double [] charge;
	double deltaW_HB[][][][];
	double expDeltaW[][][][];
	double SEGGAMMA[][];
	double SEGGAMMAPR[][][];
	
	ISegmentSolver segSolver;

	/** Flag if the Staverman-Guggenheim term is ignored */
	boolean ignoreCombinatorial = false;
	/** Flag if the residual term is ignored */
	boolean ignoreResidual = false;
	
	
	private double averageACOSMO;
	COSMOSACCompound[] comps;

	public COSMOSACMulti(){
		this(51, 2);
	}
	
	public String toString(){
		return "COSMO-MSAC";
	}
	
	/**
	 * @see #setComponents(COSMOSACCompound[])
	 */
	public COSMOSACMulti(int numberOfSegments, int numberOfDescriptors){
		this.nsegments = numberOfSegments;
		this.ndescriptors = numberOfDescriptors;
		
		this.charge = new double[nsegments];
		this.beta = new double[ndescriptors];
		this.fpol = new double[ndescriptors][ndescriptors];
		this.cHB = new double[ndescriptors][ndescriptors];
		
		double increment = (SIGMA_BOUND*2)/(nsegments-1);
		for (int i = 0; i < nsegments; i++) {
			charge[i] = -SIGMA_BOUND + increment*(double)i;
		}
		
//		// IDAC COST:0.6295141051725348, CUTOFF HB, all available IDAC
//		setResCorr(1.1495842571017678);
//		setCHB(84309.60094136275);
//		setSigmaHB(0.0071967046136033035);
//		setFpol(0.5448510600327106);
//		setIgnoreSG(false);
//		
//		// IDAC COST:0.3835742283036457, CUTOFF HB, all available IDAC without water as solute and amines
//		setResCorr(1.1280385832644837);
//		setCHB(84514.517941592);
//		setSigmaHB(0.00806375458669691);
//		setFpol(0.6123714287684722);
//		setIgnoreSG(false);
//		
//		// IDAC COST:0.5022846030794819 NP:623, CUTOFF HB, all available IDAC but amines
//		setResCorr(1.1466116315959618);
//		setCHB(37236.75882313268);
//		setSigmaHB(0.0031696732302331646);
//		setFpol(0.030186686524868365);
//		setIgnoreSG(false);
//		setAnorm(68.28177228851187);
//		setVnorm(66.69);
//		setFpol(FPOL);
		
		// IDAC COST:0.572460842572617 NP:623, CUTOFF HB, without resCorr all available IDAC but amines
//		setResCorr(1);
//		setCHB(37236.75882313268);
//		setSigmaHB(0.0031696732302331646);
//		setFpol(0.030186686524868365);
//		setIgnoreSG(false);
//		setAnorm(68.28177228851187);
//		setVnorm(66.69);

		// IDAC COST:0.5449980847647419 NP:623, CUTOFF HB, without fpol and all IDAC but amines, ketones, carbox., etc.
		setBeta(1.1166225104520564);
		setCHB(25580.016958492393);
		setSigmaHB(0.005949735966583021);
		setFpol(0.6917);
		setIgnoreSG(false);
		setAnorm(80.82815711911084);
		setVnorm(66.69);
	}

	public double getRav() {
		return rav;
	}
	public void setRav(double rav) {
		this.rav = rav;
	}
	
	public double getCoord() {
		return coord;
	}


	public void setCoord(double coord) {
		this.coord = coord;
	}
	
	/**
	 * @return the pure substance segment activity coefficient
	 */
	public double[][][] getPureSegmentGamma(){
		return SEGGAMMAPR;
	}
	/**
	 * @return the mixture segment activity coefficient
	 */
	public double[][] getMixtureSegmentGamma(){
		return SEGGAMMA;
	}
	
	/**
	 * Set the compound list for the model calculations.
	 * 
	 * <p>This function modifies the compound charge and area to match
	 * with the number of segments selected in the model.
	 * 
	 * @throws Exception 
	 * 
	 */
	public void setComponents(COSMOSACCompound []comps) throws Exception{
		this.comps = comps;
		this.ncomps = comps.length;
		
		this.VCOSMO = new double[ncomps];		

		for (int i = 0; i < comps.length; i++) {
			this.VCOSMO[i] = comps[i].Vcosmo;
			
			// Sort the compound area and charge and overwrite the compound data
//			simpleSorting(this.area[i], comps[i].charge, comps[i].area);
//			comps[i].area = this.area[i];
//			comps[i].charge = this.charge;
		}

		this.T = 300;

		z = new double[ncomps];
		for (int i = 0; i < ncomps; i++)
			z[i] = 1.0/ncomps;
		
		parametersChanged();
		setComposition(z);
	}
	
	public COSMOSACCompound[] getComps(){
		return comps;
	}
	
	/**
	 * @return the cavity volume vector
	 */
	public double[] getCavityVolume(){
		return VCOSMO;
	}

	/**
	 * Sets the system temperature in Kelvin
	 * @param T the temperature
	 */
	public void setTemperature(double T){
		this.T = T;
		this.inv_RT = 1.0 / (RGAS * T);
		
		// update this for T dependent HB
		calculeDeltaW_HB();
		
		double hbfactor = 1;
//		System.out.println("expDeltaW_RT");
		for (int m = 0; m < nsegments; m++) {
			for(int n = 0; n < nsegments; n++) {
				// classical
				double chargemn = charge[m]+charge[n];
				double chargemn2 = chargemn*chargemn;
				
				for (int d = 0; d < ndescriptors; d++) {
					
					for (int d2 = 0; d2 < ndescriptors; d2++) {
						double deltaWmn = fpol[d][d2] * (alpha/2.0)*chargemn2;
						expDeltaW[d][d2][m][n] = Math.exp(-(deltaWmn + hbfactor*deltaW_HB[d][d2][m][n]) * inv_RT);

						// System.out.println(expDeltaW[d][d2][m][n]);
					}
				}
			}
		}
		
		// ITERATION FOR SEGMENT ACITIVITY COEF (PURE SPECIES) (temperature dependent only)
		for(int i=0; i<ncomps; ++i){
			segSolver.solveMulti(comps[i].areaMulti, 1.0/ACOSMO[i], SEGGAMMAPR[i],	expDeltaW, TOLERANCE);
		}
	}
	
	public double getT(){
		return T;
	}

	/**
	 * Sets the composition of the mixture in molar basis.
	 * @param z the molar composition.
	 */
	public void setComposition(double []z){
		for (int i = 0; i < ncomps; i++)
			this.z[i] = z[i];

		// CALCULATE THE MIXTURE SIGMA PROFILE
		averageACOSMO = 0.0;
		for (int i = 0; i < ncomps; i++) {
			averageACOSMO += z[i]*ACOSMO[i];
		}
		for(int d=0; d<ndescriptors; ++d){
			for(int m=0; m<nsegments; ++m){
				double numer = 0;
				for (int i = 0; i < ncomps; i++) {
					numer += z[i]*comps[i].areaMulti[d][m]; 
				}
				PROFILE[d][m] = numer/averageACOSMO;
			}
		}
		
		// update the profiles which depend on the composition
		setTemperature(this.T);
	}
	
	public double[] getComposition(){
		return z;
	}

	/**
	 * Reallocates all mixture dependent variables.
	 * This function calculates all quantities which do not depend on composition
	 * or temperature.
	 */
	public void parametersChanged(){
		double aEff = Math.PI*rav*rav;
		aEff = 7.5;
		alpha = 0.3*Math.pow(aEff, 1.5)/E0;

		if(ACOSMO==null || ACOSMO.length!=ncomps){		
			ACOSMO = new double[ncomps];
			PROFILE = new double[ndescriptors][nsegments];
			deltaW_HB = new double[ndescriptors][ndescriptors][nsegments][nsegments];
			expDeltaW = new double[ndescriptors][ndescriptors][nsegments][nsegments];

			SEGGAMMA = new double[ndescriptors][nsegments];
			SEGGAMMAPR = new double[ncomps][ndescriptors][nsegments];
			//		segSolver = new SegmentSolverNewton();
			segSolver = new SegmentSolverSimple();

			RNORM = new double[ncomps];
			QNORM = new double[ncomps];
			li = new double[ncomps];
		}
		for(int i=0; i<ncomps; ++i){
			ACOSMO[i] = 0;
			for (int d = 0; d < ndescriptors; d++) {
				for (int m = 0; m < nsegments; m++) {
					ACOSMO[i] += comps[i].areaMulti[d][m];

					// initialize all SEGGAMMAPR (reused between calculations)
					SEGGAMMAPR[i][d][m] = 1.0;
					// initialize all SEGGAMMA (reused between calculations)
					SEGGAMMA[d][m] = 1.0;
				}
			}
		}
		// some more composition independent properties
		for(int i=0; i<ncomps; ++i){
			RNORM[i] = VCOSMO[i]/vnorm;
			QNORM[i] = ACOSMO[i]/anorm;
			li[i] = (coord/2.0)*(RNORM[i]-QNORM[i])-(RNORM[i]-1.0);
		}
		
		calculeDeltaW_HB();
	}
	
	protected void calculeDeltaW_HB(){
		for(int d=0; d<ndescriptors; ++d){
			for(int d2=0; d2<ndescriptors; ++d2){
				for(int m=0; m<nsegments; ++m){
					for(int n=0; n<nsegments; ++n){
						int ACC = n, DON = m;
						if(charge[m]>=charge[n]){
							ACC = m;
							DON = n;
						}
						// Hydrogen Bond effect:
						double chargeAcc = Math.max(0.0, charge[ACC] - sigmaHB);
						double chargeDon = Math.min(0.0, charge[DON] + sigmaHB2);
						
						// some sort of saturation
						chargeAcc = Math.min(chargeAcc, sigmaHB3);
						chargeDon = Math.max(chargeDon, -sigmaHB3);
						
						double hb;
						switch (hbType) {
						case HBTYPE_HSIEH:
							hb = chargeAcc-chargeDon;
							hb *= hb;
							break;
						case HBTYPE_FIXED:
							// is always the same (just a scaling factor to get the same order of magnitude)
							hb = (0.010*0.010);
							break;
						case HBTYPE_COVALENT:
//							hb = 0.006*chargeAcc;
//							chargeDon = 0.010 + chargeDon;
//							hb = 0.006*chargeDon;
							hb = chargeAcc-chargeDon;
							hb *= hb;
							if(chargeDon < -0.007 || chargeAcc > 0.010)
								hb = 0;
							break;
						default:
							hb = chargeAcc*chargeDon;
							break;
						}
						
						// Klamt, Fluid Phase Equilib. 2000
						// double cHBT_c = 1.5;
						double cHBT = 1; // Math.max(0, 1 + cHBT_c * (298.15/T - 1));
						
						hb = -Math.abs(hb);

						deltaW_HB[d][d2][m][n] = cHB[d][d2]*cHBT* hb;
					}
				}
			}
		}
	}

	/**
	 * @return the current segment solver.
	 */
	public ISegmentSolver getSegSolver() {
		return segSolver;
	}

	/**
	 * Set the segment solver, default is {@link SegmentSolverNewton}.
	 * @param segSolver the segment solver to be used
	 */
	public void setSegSolver(ISegmentSolver segSolver) {
		this.segSolver = segSolver;
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
		// Solve for the SEGGAMMA of the mixture
		segSolver.solveMulti(PROFILE, 1.0, SEGGAMMA, expDeltaW, TOLERANCE);
		
		// THE STAVERMAN-GUGGENHEIM EQUATION
		double sum_zi_qi = 0;
		double sum_zi_ri = 0;
		double sum_ziR_ri = 0;
		for(int i=0; i<ncomps; ++i){
			sum_zi_qi += z[i]*QNORM[i];
			sum_zi_ri += z[i]*RNORM[i];
			sum_ziR_ri += z[i]*Math.pow(RNORM[i], rPower);
		}
		
		double aEff = Math.PI*rav*rav;
		aEff = 7.5;

		double sum_zi_li = 0;
		for(int i=0; i<ncomps; ++i)
			sum_zi_li += z[i]*li[i];

		for(int i=0; i<ncomps; ++i){
			double phi_zi = RNORM[i]/sum_zi_ri;
			double psi_zi = Math.pow(RNORM[i], rPower)/sum_ziR_ri;
			double theta_zi = QNORM[i]/sum_zi_qi;

			double lnGammaSG = 0;
			if(!ignoreCombinatorial){
				// original expression
//				double theta_phi = theta_zi/phi_zi;
//				lnGammaSG += Math.log(phi_zi)+
//				(coord/2)*QNORM[i]*Math.log(theta_phi)
//				+ li[i] - (phi_zi)*sum_zi_li;

				// expression from DOI:10.1021/ie070465z
				double phi_theta = phi_zi/theta_zi;
				lnGammaSG += Math.log(psi_zi) + 1 - psi_zi - (coord/2)*QNORM[i]*
					( Math.log(phi_theta) + 1 - phi_theta);

			}
			
			// CALCULATION OF GAMMAS
			double lnGammaRestoration = 0.0;
			if(!ignoreResidual){
				for (int d = 0; d < ndescriptors; d++) {
					for(int m=0; m<nsegments; ++m){
						// avoid computing for null areas
						if(comps[i].areaMulti[d][m]==0)
							continue;

						// lnGammaRestoration += (sigma[i][m]/aEffPrime)*(Math.log(SEGGAMMA[m]/(SEGGAMMAPR[i][m])));
						double lnMixSeg = Math.log(SEGGAMMA[d][m]);
						double lnMixPure = Math.log(SEGGAMMAPR[i][d][m]);

						lnGammaRestoration += beta[d]*(comps[i].areaMulti[d][m]/aEff)*(lnMixSeg - lnMixPure);
					}
				}
			}
			lnGama[i+start] = lnGammaRestoration + lnGammaSG;
		}
	}

	public boolean isIgnoreSG() {
		return ignoreCombinatorial;
	}
	public boolean isIgnoreResidual() {
		return ignoreResidual;
	}

	/**
	 * Sets if the Staverman-Guggenheim term should be ignored.
	 * @param ignoreSG <code>true</code> if the term should be ignored
	 */
	public void setIgnoreSG(boolean ignoreSG) {
		this.ignoreCombinatorial = ignoreSG;
	}
	
	/**
	 * Sets if the residual term should be ignored.
	 * @param ignoreResidual <code>true</code> if the term should be ignored
	 */
	public void setIgnoreResidual(boolean ignoreResidual) {
		this.ignoreResidual = ignoreResidual;
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

	public double getRPower() {
		return rPower;
	}
	public void setRPower(double rPower) {
		this.rPower = rPower;
	}


	public double getCHB() {
		return cHB[0][0];
	}
	public double getCHB(int i, int j) {
		return cHB[i][j];
	}


	public void setCHB(double chb) {
		for (int i = 0; i < ndescriptors; i++) {
			for (int j = 0; j < ndescriptors; j++) {
				this.cHB[i][j] = Math.abs(chb);
			}
		}
	}
	public void setCHB(int i, int j, double chb) {
		this.cHB[i][j] = Math.abs(chb);
		this.cHB[j][i] = Math.abs(chb);
	}


	public double getSigmaHB() {
		return sigmaHB;
	}
	public void setSigmaHB(double sigmaHB) {
		this.sigmaHB = Math.abs(sigmaHB);
	}
	public double getSigmaHB2() {
		return sigmaHB2;
	}
	public void setSigmaHB2(double sigmaHB) {
		this.sigmaHB2 = sigmaHB;
	}
	public double getSigmaHB3() {
		return sigmaHB3;
	}
	public void setSigmaHB3(double sigmaHB) {
		this.sigmaHB3 = sigmaHB;
	}

	
	public double getFpol() {
		return fpol[0][0];
	}
	public double getFpol(int i, int j) {
		return fpol[i][j];
	}
	public void setFpol(double fpol) {
		for (int i = 0; i < ndescriptors; i++) {
			for (int j = 0; j < ndescriptors; j++) {
				this.fpol[i][j] = Math.abs(fpol);				
			}
		}
	}
	
	/**
	 * Set the fpol value for a given descriptor
	 * @param i the descriptor
	 * @param fpol the value
	 */
	public void setFpol(int i, int j, double fpol) {
		this.fpol[i][j] = Math.abs(fpol);
		this.fpol[j][i] = Math.abs(fpol);
	}
	
	public double getBeta() {
		return beta[0];
	}
	public double getBeta(int i) {
		return beta[i];
	}	
	/**
	 * Set the beta value for all descriptors
	 * @param beta the value
	 */
	public void setBeta(double beta) {
		for (int i = 0; i < this.beta.length; i++) {
			this.beta[i] = beta;
		}
	}
	
	/**
	 * Set the beta parameter for a given descriptor
	 * @param i the descriptor
	 * @param beta the value
	 */
	public void setBeta(int i, double beta) {
		this.beta[i] = beta;
	}

	public int getHbType() {
		return hbType;
	}

	public void setHbType(int hbType) {
		this.hbType = hbType;
	}
}
