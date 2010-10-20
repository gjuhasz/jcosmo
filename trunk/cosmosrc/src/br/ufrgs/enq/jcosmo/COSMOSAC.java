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
	double rPower = 1;

	/** Default coordination number (KLAMT USED 7.2) */
	public static final double COORD = 10.0;
	double coord = COORD;

	public static final double SIGMAHB = 0.0084;
	public static final double CHB = 85580.0;
	
	double sigmaHB = SIGMAHB;
	double sigmaHB2 = SIGMAHB;
	double sigmaHB3 = 0;
	double cHB = CHB;
	
	double sigmaDisp = 0.002;
	double cDisp = 1000;

	private static double SIGMA_BOUND = 0.025;
	
	double alpha;

	private static final double TOLERANCE = 1e-8;

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
	
	double permitPure[];
	double permit;
	
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
		
//		// article results
//		setBeta(1.12);
//		setCHB(25580.016958492393);
//		setSigmaHB(0.00595);
//		setSigmaHB2(getSigmaHB());
//		setSigmaHB3(1.0);
//		setFpol(FPOL); // ?
//		setIgnoreSG(false);
//		setCoord(10.0);
//		setAnorm(80.83);
//		setVnorm(66.69);
		
		
		setAnorm(124.18305909611365);
		setRPower(0.6372590453541258);
		setBeta(1.399460298613211);
		setFpol(0.73377065315896);
		setCHB(0);
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
		this.ACOSMO = new double[ncomps];
		this.permitPure = new double[ncomps];
		this.area = new double[ncomps][nsegments];

		for (int i = 0; i < comps.length; i++) {
			if(comps[i].area == null)
				comps[i] = COSMOSACDataBase.getInstance().getComp(comps[i].name);
			this.VCOSMO[i] = comps[i].Vcosmo;
			this.area[i] = comps[i].area;
			
			ACOSMO[i] = 0;
			double s2 = 0;
			double []charge = comps[i].charge;
			for (int m = 0; m < charge.length; m++) {
				ACOSMO[i] += area[i][m];
				s2 += area[i][m]*charge[m]*charge[m];
			}
			s2 /= ACOSMO[i];
			s2*=10000;
			// from linear fit of some substances, with a lower limit
			permitPure[i] = Math.max(1.8, s2*71.3-7.74);
			
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
		
//		System.out.println("expDeltaW_RT");
		for (int m = 0; m < nsegments; m++) {
			for(int n = 0; n < nsegments; n++) {
				double scale = 1;

				// pure expDeltaW
				for (int i = 0; i < ncomps; i++) {
//					scale = (permitPure[i]-1)/(permitPure[i]+0.5)/0.98;
//					scale = Math.cbrt(scale);
					expDeltaWPure[i][m][n] = Math.exp(-(deltaW[m][n] + deltaW_HB[m][n])*scale * inv_RT);
				}
				
				// mixture expDeltaW
//				scale = ((permit-1)/(permit+0.5))/0.98;
//				scale = Math.cbrt(scale);
				expDeltaW[m][n] = Math.exp(-(deltaW[m][n] + deltaW_HB[m][n])*scale * inv_RT);
				
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
		
		double s2 = 0;
		double []charge = comps[0].charge;
		for (int m = 0; m < PROFILE.length; m++) {
			s2 += PROFILE[m]*charge[m]*charge[m];
		}
		s2*=10000;
		// from linear fit of some substances, with a lower limit
		permit = Math.max(1.8, s2*71.3-7.74);
		
//		permit = 0;
//		for (int i = 0; i < ncomps; i++) {
//			permit += 1/(z[i]*ACOSMO[i]*permitPure[i]);
//		}
//		permit = 1/(permit*averageACOSMO);
		
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
			for (int m = 0; m < area[i].length; m++) {
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
				// original formulation
				double chargemn = charge[m]+charge[n];
				deltaW[m][n] = (fpol*alpha/2.0)*chargemn*chargemn;
				
//				// only positive charges are polarized
//				double chargem = charge[m];
//				double chargen = charge[n];
//				if(chargem < 0)
//					chargem *= fpol;
//				if(chargen < 0)
//					chargen *= fpol;
//				double chargemn = chargem+chargen;
//				deltaW[m][n] = (alpha/2.0)*chargemn*chargemn;
			}
		}
	}
	
	protected void calculeDeltaW_HB(){
		for(int m=0; m<nsegments; ++m){

			for(int n=0; n<nsegments; ++n){
				double hb = 0;
				double chargeacc = charge[m];
				double chargedon = charge[n];

				// Hydrogen Bond effect:
				if(chargeacc < chargedon){
					double tmp = chargedon;
					chargedon = chargeacc;
					chargeacc = tmp;
				}
				
				// a saturation on this effect
//				chargedon = Math.max(-sigmaHB3, chargedon);
//				chargeacc = Math.min( sigmaHB3, chargeacc);
				
				hb = Math.max(0.0, chargeacc - sigmaHB)*Math.min(0.0, chargedon + sigmaHB2);
//				hb = Math.min(0.0, charge[ACC]*charge[DON] - sigmaHB*sigmaHB);
				hb = -Math.abs(hb);
				
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
			if(!ignoreSG){
				// original expression
//				double theta_phi = theta_zi/phi_zi;
//				lnGammaSG += Math.log(phi_zi)+
//				(coord/2)*QNORM[i]*Math.log(theta_phi)
//				+ li[i] - (phi_zi)*sum_zi_li;

				// expression from DOI:10.1021/ie070465z
				double phi_theta = phi_zi/theta_zi;
				lnGammaSG += Math.log(psi_zi) + 1 - psi_zi - (coord/2)*QNORM[i]*
					( Math.log(phi_theta) + 1 - phi_theta);
				
				// combinatorial from Klamt, 1998 (eq 30. DOI 10.1021/jp980017s)
				// there is a major problem with this expression (to be published)
				// lnGammaSG = 0.14*Math.log(sum_zi_qi/QNORM[i]);
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

	public double getRPower() {
		return rPower;
	}

	public void setRPower(double rPower) {
		this.rPower = rPower;
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

	public void setCDisp(double cDisp) {
		this.cDisp = cDisp;
	}
	public double getCDisp() {
		return cDisp;
	}
	
	public double getSigmaDisp() {
		return sigmaDisp;
	}
	public void setSigmaDisp(double sigmaDisp) {
		this.sigmaDisp = sigmaDisp;
	}


	public double getSigmaHB() {
		return sigmaHB;
	}
	public void setSigmaHB(double sigmaHB) {
		this.sigmaHB = sigmaHB;
		this.sigmaHB2 = sigmaHB;
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
