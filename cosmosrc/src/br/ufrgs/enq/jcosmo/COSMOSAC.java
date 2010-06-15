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
 * <li>It can evaluate infinite dilution (x=0) activity coefficients (original code
 * had a bug and returned NaN for this case)
 * <li>There are some changes on the segment activity converging mechanism (much better performance)
 * </ul>
 * 
 * @author Rafael de Pelegrini Soares
 * 
 */
public class COSMOSAC {

	static final double E0 = 2.395e-4;
	public static final double FPOL = 0.6917; // Ind. Eng. Chem. Res., Vol. 46, No. 22, 2007
	/** Default averaging radius */
	public static final double RAV = 0.81764200000000;

	double rav = RAV;
	double fpol = FPOL;
	double beta = 1;

	static final double RGAS = 0.001987; // kcal/mol/K
	public static final double VNORM = 66.69;
	public static final double ANORM = 79.53;

	double vnorm = VNORM;
	double anorm = ANORM;

	/** Default coordination number (KLAMT USED 7.2) */
	public static final double COORD = 10.0;
	double coord = COORD;

	public static final double SIGMAHB = 0.0084;
	public static final double CHB = 85580.0;
	
	double sigmaHB = SIGMAHB;
	double sigmaHB2 = 0;
	double sigmaHB3 = 1;
	double cHB = CHB;

	private static double SIGMA_BOUND = 0.025;
	
	double alpha;

	private static final double TOLERANCE = 1e-6;

	protected double T;
	protected double inv_RT;
	protected double[] z;

	/** The number of compounds in the mixture */
	int ncomps;
	/** The number of independent segments to be considered by the model */ 
	int nsegments;

	/// Cavity volume for each substance
	double [] VCOSMO;
	double [] ACOSMO;

	double [] RNORM;
	double [] QNORM;
	double [] li;

	double [] PROFILE;

	double [][] area;
	double [] charge;
	double deltaW[][];
	double deltaW_HB[][];
	double expDeltaWPure[][][];
	double expDeltaW[][];
	double SEGGAMMA[];
	double SEGGAMMAPR[][];
	
	ISegmentSolver segSolver;

	/** Flag if the Staverman-Guggenheim term is ignored */
	boolean ignoreSG = false;
	private double averageACOSMO;
	COSMOSACCompound[] comps;

	public COSMOSAC(){
		this(51);
	}
	
	public String toString(){
		return "COSMO-SAC(DMol3)";
	}
	
	/**
	 * @see #setComponents(COSMOSACCompound[])
	 */
	COSMOSAC(int numberOfSegments){
		this.nsegments = numberOfSegments;
		this.charge = new double[nsegments];
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
		
		setBeta(1.0947446969598862);
		setCHB(35641.7258589243);
		setSigmaHB(0.002621822743558675);
		setSigmaHB2(0.0);
		setSigmaHB3(1.0);
		setFpol(0.024526693237433943);
		setIgnoreSG(false);
		setCoord(10.0);
		setAnorm(71.72908857440859);
		setVnorm(66.69);
		
		
		// article results
		setBeta(1.12);
		setCHB(25580.016958492393);
		setSigmaHB(0.00595);
		setSigmaHB2(0.0);
		setSigmaHB3(1.0);
		setFpol(FPOL); // ?
		setIgnoreSG(false);
		setCoord(10.0);
		setAnorm(80.83);
		setVnorm(66.69);
		
		// electrostatic deltaW, COST:0.5122766186327498 NP:672
		setBeta(1.1983740582107176);
		setCHB(35733.56320058883);
		setSigmaHB(0.0029985190152935796);
		setSigmaHB2(1.0);
		setSigmaHB3(1.0);
		setFpol(0.07497944979786468);
		setIgnoreSG(false);
		setCoord(10.0);
		setAnorm(80.83);
		setVnorm(66.69);
		
//		// idac/nonaqueous.csv AARD:0.36629188907623855 NP:361
//		// idac/aqueous.csv AARD:0.7007706720246113 NP:311
//		// COST:0.5210874756068953
//		setBeta(1.1476530578891553);
//		setCHB(31054.041024247716);
//		setSigmaHB(0.0033043145713630077);
//		setSigmaHB2(0.48530335590723794);
//		setSigmaHB3(1.0);
//		setFpol(0.45471090287246685);
//		setIgnoreSG(false);
//		setCoord(10.0);
//		setAnorm(64.7412501912964);
//		setVnorm(66.69);

		// solvent water only COST:0.4510264290028932
		setBeta(1.2024778569098764);
		setCHB(43806.493230652006);
		setSigmaHB(0.003518327404610549);
		setSigmaHB2(0.0);
		setSigmaHB3(1.0);
		setFpol(0.07742269047560099);
		setIgnoreSG(false);
		setCoord(10.0);
		setAnorm(80.21650597414602);
		setVnorm(66.69);
		
		setCHB(43806.493230652006);
		
//		setSigmaHB(0.001);
//		setSigmaHB3(0.001);
//		setCHB(22000);
//		setFpol(0.1);
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
	public double[][] getPureSegmentGamma(){
		return SEGGAMMAPR;
	}
	/**
	 * @return the mixture segment activity coefficient
	 */
	public double[] getMixtureSegmentGamma(){
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
		this.area = new double[ncomps][nsegments];

		for (int i = 0; i < comps.length; i++) {
			this.VCOSMO[i] = comps[i].Vcosmo;
			this.area[i] = comps[i].area;
			
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
	
	/**
	 * Sorts a given sigma profile to match the model charge discretization.
	 * 
	 * @param sortedArea will contain the resulting sorted area
	 * @param unsortedSigma the compound unsorted charge list
	 * @param unsortedArea the unsorted area of each charge
	 */
	void simpleSorting(double sortedArea[], double unsortedSigma[], double[]unsortedArea){
		double increment = (SIGMA_BOUND*2)/(double)(nsegments-1);
		double unsortedIncrement = unsortedSigma[1]-unsortedSigma[0];

		// SIGMA PROFILE SORTING TAKEN FROM LIN DISSERTATION**
		for (int J = 0; J < unsortedSigma.length; J++) {
			// out of the bounds
			if(unsortedSigma[J]<=-SIGMA_BOUND){
				sortedArea[0] += unsortedArea[J]*unsortedIncrement/increment;
				continue;
			}
			if(unsortedSigma[J]>=SIGMA_BOUND){
				sortedArea[sortedArea.length-1] += unsortedArea[J]*unsortedIncrement/increment;
				continue;
			}
			
			int pos = (int)((unsortedSigma[J]-charge[0])/increment);
			// Each point represents the center of an interval, so we distribute it
			// in the two adjacent points
			double unsortedStart = unsortedSigma[J]-unsortedIncrement/2;
			double unsortedEnd = unsortedSigma[J]+unsortedIncrement/2;
			double posEnd = charge[pos]+increment/2;
			
			sortedArea[pos]+= unsortedArea[J]/increment * Math.max(0, posEnd - unsortedStart);
			sortedArea[pos+1]+= unsortedArea[J]/increment * Math.max(0, unsortedEnd - posEnd);
		}
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

				// pure expDeltaW
				for (int i = 0; i < ncomps; i++) {
					expDeltaWPure[i][m][n] = Math.exp(-(deltaW[m][n] + hbfactor*deltaW_HB[m][n]) * inv_RT);
				}
				
				// mixture expDeltaW
				expDeltaW[m][n] = Math.exp(-(deltaW[m][n] + hbfactor*deltaW_HB[m][n]) * inv_RT);
				
//				System.out.println(expDeltaW[m][n]);
			}
		}
		
		// ITERATION FOR SEGMENT ACITIVITY COEF (PURE SPECIES) (temperature dependent only)
		for(int i=0; i<ncomps; ++i){
			segSolver.solve(area[i], 1.0/ACOSMO[i], SEGGAMMAPR[i],	expDeltaWPure[i], TOLERANCE);
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
		for(int m=0; m<nsegments; ++m){
			double numer = 0;
			for (int i = 0; i < ncomps; i++) {
				numer += z[i]*area[i][m]; 
			}
			PROFILE[m] = numer/averageACOSMO;
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
		
		ACOSMO = new double[ncomps];
		PROFILE = new double[nsegments];
		deltaW = new double[nsegments][nsegments];
		deltaW_HB = new double[nsegments][nsegments];
		expDeltaW = new double[nsegments][nsegments];
		expDeltaWPure = new double[ncomps][nsegments][nsegments];
		
		SEGGAMMA = new double[nsegments];
		SEGGAMMAPR = new double[ncomps][nsegments];
//		segSolver = new SegmentSolverNewton();
		segSolver = new SegmentSolverSimple();
		
		RNORM = new double[ncomps];
		QNORM = new double[ncomps];
		li = new double[ncomps];
		for(int i=0; i<ncomps; ++i){
			ACOSMO[i] = 0;
			for (int m = 0; m < area[i].length; m++) {
				ACOSMO[i] += area[i][m];

				// initialize all SEGGAMMAPR (reused between calculations)
				SEGGAMMAPR[i][m] = 1.0;
			}
		}
		// some more composition independent properties
		for(int i=0; i<ncomps; ++i){
			RNORM[i] = VCOSMO[i]/vnorm;
			QNORM[i] = ACOSMO[i]/anorm;
			li[i] = (coord/2.0)*(RNORM[i]-QNORM[i])-(RNORM[i]-1.0);
		}
		
		calculeDeltaW();
		calculeDeltaW_HB();
	}
	
	protected void calculeDeltaW(){
		for(int m=0; m<nsegments; ++m){
			// initialize all SEGGAMMA (reused between calculations)
			SEGGAMMA[m] = 1.0;

			for(int n=0; n<nsegments; ++n){
//				// electrostatic formulation
//				double chargemn = charge[m]*charge[n];				
//				// with a minimum cutt-off
////				chargemn = Math.max(-0.018, chargemn);
////				chargemn = Math.min(0, chargemn);
//				if(chargemn>0)
//					chargemn*=1.4;
//				else
//					chargemn/=1.4;
//				deltaW[m][n] = (fpol*alpha/2.0)*chargemn;
				
				// original formulation
//				double chargemn = charge[m]+charge[n];
//				deltaW[m][n] = (fpol*alpha/2.0)*chargemn*chargemn;
				
				double chargemn = charge[m]*charge[n];
				if(chargemn<0)
					deltaW[m][n] = (fpol*alpha/2.0)*chargemn;
				else
					deltaW[m][n] = (fpol*alpha/2.0)*chargemn*sigmaHB2;

//				double chargem = Math.abs(charge[m]);
//				double chargen = Math.abs(charge[n]);
//				double chargemn = chargem+chargen;
//				if(chargem<sigmaHB2 || chargen<sigmaHB2){
//					// attraction between non-polar segments
//					deltaW[m][n] = -(fpol*alpha/2.0)*chargemn*chargemn;
//					
////					chargem+=sigmaHB3;
////					chargen+=sigmaHB3;
////					deltaW[m][n] = -(fpol*alpha/2.0)*(chargem*chargem + chargen*chargen);
//					
//					// attraction between non-polar segments
////					double inductor = Math.max(chargem,	chargen);
////					deltaW[m][n] = -(fpol*alpha/2.0)*inductor;
//				}
////				else{
////					chargemn = charge[m]+charge[n];
////					deltaW[m][n] = (fpol*alpha/2.0)*chargemn*chargemn;
////				}
			}
		}
	}
	
	protected void calculeDeltaW_HB(){
		for(int m=0; m<nsegments; ++m){

			for(int n=0; n<nsegments; ++n){
				double hb = 0;

				// Hydrogen Bond effect:
				int ACC = n, DON = m;
				if(charge[m]>=charge[n]){
					ACC = m;
					DON = n;
				}
				hb = Math.max(0.0, charge[ACC] - sigmaHB)*Math.min(0.0, charge[DON] + sigmaHB);
//				hb = Math.min(0.0, charge[ACC]*charge[DON] - sigmaHB*sigmaHB);
				hb = -Math.abs(hb);
				
//				 simple electrostatics
//				if(Math.abs(charge[m])>sigmaHB && Math.abs(charge[n])>sigmaHB)
//					hb = charge[m]*charge[n];
				
				// Klamt, Fluid Phase Equilib. 2000
//				double cHBT_c = 1.5;
				double cHBT = 1; // Math.max(0, 1 + cHBT_c * (298.15/T - 1));
				
				deltaW_HB[m][n] = cHB*cHBT* hb;
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
		segSolver.solve(PROFILE, 1.0, SEGGAMMA, expDeltaW, TOLERANCE);
		
		// THE STAVERMAN-GUGGENHEIM EQUATION
		double BOTTHETA = 0;
		double BOTPHI = 0;
		for(int i=0; i<ncomps; ++i){
			BOTTHETA += z[i]*QNORM[i];
			BOTPHI += z[i]*RNORM[i];
		}
		
		double aEff = Math.PI*rav*rav;
		aEff = 7.5;

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
			for(int m=0; m<nsegments; ++m){
//				lnGammaRestoration += (sigma[i][m]/aEffPrime)*(Math.log(SEGGAMMA[m]/(SEGGAMMAPR[i][m])));
				double lnMixSeg = Math.log(SEGGAMMA[m]);
				double lnMixPure = Math.log(SEGGAMMAPR[i][m]);
				
				lnGammaRestoration += (area[i][m]/aEff)*(lnMixSeg - lnMixPure);
			}
			lnGama[i+start] = beta*lnGammaRestoration + lnGammaSG;
		}
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
		return fpol;
	}
	public void setFpol(double fpol) {
		this.fpol = fpol;
	}
	
	public double getBeta() {
		return beta;
	}
	public void setBeta(double beta) {
		this.beta = beta;
	}
}
