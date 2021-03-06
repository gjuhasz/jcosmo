/*
 * Copyright 2008-2010 Universidade Federal do Rio Grande do Sul
 * 
 * License should still be decided.
 * 
 * Author: Rafael de Pelegrini Soares (rafael@enq.ufrgs.br)
 */

package br.ufrgs.enq.jcosmo2;

import java.io.File;

import org.netlib.blas.BLAS;

import br.ufrgs.enq.jcosmo.SigmaProfileGenerator;
import br.ufrgs.enq.jcosmo.SigmaProfileGenerator.FileType;




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
public class COSMOSAC2 {

	static final double E0 = 2.395e-4;
	public static final double FPOL = 0.6917; // Ind. Eng. Chem. Res., Vol. 46, No. 22, 2007
	/** Default averaging radius */
	public static final double RAV = 0.81764200000000;
	
	// SIGMA compression stuff
	private static final int COMPRESSED_SEGMENTS = 61;
	private static final double SIGMA_MAX = 0.030;
	private static final double SIGMA_THRESHOLD = 2*SIGMA_MAX/(COMPRESSED_SEGMENTS-1);
	private static final double[] SIGMA_COMPRESSED = new double[COMPRESSED_SEGMENTS];
	private static final double[] AREA_COMPRESSED = new double[COMPRESSED_SEGMENTS];
	private static final double SIGMA_INCREMENT = (SIGMA_MAX*2)/(COMPRESSED_SEGMENTS-1);
	static{
		for (int m = 0; m < COMPRESSED_SEGMENTS; m++) {
			SIGMA_COMPRESSED[m] = -SIGMA_MAX + SIGMA_INCREMENT*(double)m;
		}
	}

	double fpol = FPOL;
	double beta = 1;
	
	double rav = RAV;
	double rav2 = 1.5*RAV;
	double f_ortho = 0.79209;
	double fcorr = -2.0;

	static final double RGAS = 0.001987; // kcal/mol/K
	public static final double VNORM = 66.69;
	public static final double ANORM = 79.53;

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
	
	double alpha;

	private static final double TOLERANCE = 1e-8;

	protected double T;
	protected double inv_RT;
	protected double[] z;

	/** The number of compounds in the mixture */
	int ncomps;

	double [] RNORM;
	double [] QNORM;

	double [][] PROFILE;

	double SEGGAMMA[][];
	double SEGGAMMAPR[][];
	
	/** Flag if the Staverman-Guggenheim term is ignored */
	boolean ignoreSG = false;
	private double averageACOSMO;
	
	boolean skipAverage = false;
	
	Compound[] comps;
	
	String folder = "profiles/RM1_1.18/";
	String extension = ".cos";
	FileType type = SigmaProfileGenerator.FileType.MOPAC;

	public String toString(){
		return "COSMO-SAC2";
	}
	
	/**
	 * @see #setComponents(Compound[])
	 */
	COSMOSAC2(){
		rav = 1;
		rav2 = 1.5*rav;
		f_ortho = 0.79209;
		fcorr = -2.0;

		setBeta(2.03472);
		setFpol(0.685231);
		setAnorm(64.61423251936094);
		setRPower(0.7700935807819056);
		setCHB(0);

		folder = "profiles/RM1_1.18/";
		skipAverage = false;
		rav = 1;
		rav2 = 1.5*rav;
		f_ortho = 0.8033901717857764;
		fcorr = -2.0;
		setBeta(1.988614770236198);
		setFpol(0.6776862764908341);
		
		folder = "profiles/RM1_1.18/";
		skipAverage = false;
		rav = 1.5;
		rav2 = 1.5*rav;
		f_ortho = 0.6665617963244067;
		setBeta(1.7897260990915678);
		setFpol(0.8618302987294044);
		setFcorr(-0.03928437537005114);
		
		// nonHB COST:0.1149666078072233 NP:237
		folder = "profiles/RM1_1.18/";
		skipAverage = false;
		rav = 1.5;
		rav2 = 2*rav;
		f_ortho = 0.41830198448203804;
		fcorr = -2.0;
		setBeta(1.3674741777731985);
		setFpol(1.0982211235904975);
		setFcorr(-2.859963039753584);
		
		// nonHB COST:0.14219163689050052 NP:272		
		folder = "profiles/RM1/";
		skipAverage = false;
		rav = 1.5;
		rav2 = 2*rav;
		f_ortho = 0.6057217795048631;
		setBeta(2.4226621419936434);
		setFpol(0.38473549097776416);
		setFcorr(1);
		
		folder = "profiles/RM1_1.18/";
		// linear relation from tests
		rav = 1.1;
		// eq for av
		f_ortho = rav*-0.312370 + 1.05062; // av
//		f_ortho = rav*-0.346096 + 1.126528; // av3

		
//		// nonHB COST:0.1404634741935444
//		folder = "profiles/RM1/";
//		skipAverage = true;
//		setBeta(1.668791603688795);
//		setFpol(0.5527341346636111);
	}

	public double getFcorr() {
		return fcorr;
	}

	public void setFcorr(double fcorr) {
		this.fcorr = fcorr;
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
	public double[][] getMixtureSegmentGamma(){
		return SEGGAMMA;
	}
	
	protected void loadCompoundData(Compound comp) throws Exception {
		
		String name = comp.name.replace(' ','_').toUpperCase();
		File file = new File(folder + name + extension);
		if(!file.exists())
			name = name.replace('-','_');
			
		SigmaProfileGenerator s = new SigmaProfileGenerator(type , this.rav, 0);
		s.parseFile(folder + name + extension);
		comp.area = s.getOriginalArea();
		comp.x = s.getX();
		comp.y = s.getY();
		comp.z = s.getZ();
		
		comp.volumeTotal = s.getVolume();
		comp.areaTotal = 0;
		
		if(skipAverage){
			// Either use the original charge density (no averaging)
			comp.sigmaAvg = s.getOriginalChargeDensity();
			for (int m = 0; m < comp.area.length; m++)
				comp.areaTotal += comp.area[m];
		}
		else{
			// Or use the averaged charges
			comp.sigmaAvg = s.getAveragedChargeDensity();
			//			SigmaProfileGenerator s2 = new SigmaProfileGenerator(type, this.rav2, 0);
			//			s2.parseFile(folder + name + extension);
			//			comp.sigmaOrtho = s2.getAveragedChargeDensity();
			comp.sigmaDelta = s.getOriginalChargeDensity();

			double energy0 = 0, energyAvg = 0;
			for (int m = 0; m < comp.area.length; m++){
				energy0 += comp.area[m] * comp.sigmaDelta[m]*comp.sigmaDelta[m];
				energyAvg += comp.area[m] * comp.sigmaAvg[m]*comp.sigmaAvg[m];
				
				comp.areaTotal += comp.area[m];
				comp.sigmaDelta[m] = f_ortho*comp.sigmaDelta[m] - comp.sigmaAvg[m];
			}
//			double [] sigmaDeltaAvg = new double[comp.sigmaDelta.length];
//			SigmaProfileGenerator2.averageCharges(rav, comp.sigmaDelta, sigmaDeltaAvg, comp.area, comp.x, comp.y, comp.z);
//			comp.sigmaDelta = sigmaDeltaAvg;
			
			double sigmaScale = Math.sqrt(energy0/energyAvg);
			for (int m = 0; m < comp.area.length; m++){
				comp.sigmaAvg[m] *= sigmaScale;
				comp.sigmaDelta[m] = 0;
			}
		}
		
		compressProfile(comp);
	}
	
	/**
	 * Compress the sigma profile of the given compound.
	 * 
	 * <p>Only nonzero elements are kept and those with only a small
	 * deviation correction are aggregated as possible.
	 * 
	 * @param comp
	 */
	private void compressProfile(Compound comp){
		int last = AREA_COMPRESSED.length-1;
		int nremain = comp.area.length, nsorted = 0;
		
		for (int J = 0; J < comp.area.length; J++) {
			if(comp.sigmaDelta!=null && Math.abs(comp.sigmaDelta[J]) > SIGMA_THRESHOLD)
				continue;
			
			int TMP = (int)((comp.sigmaAvg[J]+SIGMA_MAX)/SIGMA_INCREMENT);
			
			if(TMP<0 || TMP > last)
				continue;
			
			// Charges are simply put into a range 
			if(comp.sigmaAvg[J]-SIGMA_COMPRESSED[TMP] > SIGMA_INCREMENT/2)
				TMP++;

			if(AREA_COMPRESSED[TMP]==0)
				nsorted++;
			AREA_COMPRESSED[TMP] += comp.area[J];
			comp.area[J] = 0;
			nremain--;
		}
		
		// assemble the new vectors
		int newLength = nsorted + nremain;
		double newArea[] = new double[newLength];
		double newSigma[] = new double[newLength];
		double newSigmaOrtho[] = new double[newLength];
		int index = 0;
		// add the sorted elements
		for (int i = 0; i < AREA_COMPRESSED.length; i++) {
			if(AREA_COMPRESSED[i] == 0)
				continue;
			newArea[index] = AREA_COMPRESSED[i];
			AREA_COMPRESSED[i] = 0; //reset for the next use
			newSigma[index] = SIGMA_COMPRESSED[i];
			index++;
		}
		// add the remaining elements
		for (int i = 0; i < comp.area.length; i++) {
			if(comp.area[i] == 0)
				continue;
			newArea[index] = comp.area[i];
			newSigma[index] = comp.sigmaAvg[i];
			newSigmaOrtho[index] = comp.sigmaDelta[i];
			index++;
		}
		
		// replace the old vectors
		comp.area = newArea;
		comp.sigmaAvg = newSigma;
		comp.sigmaDelta = newSigmaOrtho;
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
	public void setComponents(Compound []comps) throws Exception{
		this.comps = comps;
		this.ncomps = comps.length;
		
		
		this.PROFILE = new double[ncomps][];
		this.SEGGAMMA = new double[ncomps][];
		this.SEGGAMMAPR = new double[ncomps][];

		for (int i = 0; i < comps.length; i++) {
			loadCompoundData(comps[i]);
			
			int length = comps[i].area.length;
			this.PROFILE[i] = new double[length];
			this.SEGGAMMA[i] = new double[length];
			this.SEGGAMMAPR[i] = new double[length];
		}
		
		z = new double[ncomps];
		for (int i = 0; i < ncomps; i++)
			z[i] = 1.0/ncomps;
		
		parametersChanged();
		setComposition(z);
		setTemperature(300);
	}
	
	public Compound[] getComps(){
		return comps;
	}
	
	/**
	 * Sets the system temperature in Kelvin
	 * @param T the temperature
	 */
	public void setTemperature(double T){
		this.T = T;
		this.inv_RT = 1.0 / (RGAS * T);
		
		// with the temperature we can now compute the segment activity of the pure substances:
		solveGamma(true, SEGGAMMAPR, TOLERANCE);
	}
	
	/**
	 * Solve the iterative equations for the segment activity coefficients.
	 * 
	 * @param pure if calculation should consider substances as pure
	 * @param segGamma the resulting segment activity coefficients
	 * @param tol the tolerance to be considered
	 */
	private void solveGamma(boolean pure, double[][] segGamma, double tol) {
		
		int niter = 0;
		int maxiter = 100;
		double norm = -1;
		BLAS blas = BLAS.getInstance();
		
		while(true){
			for (int i = 0; i < ncomps; i++) {
				Compound compi = comps[i];
				double []segGammai = segGamma[i];
				
				double factor = pure ? comps[i].areaTotal : 1.0;
				
				for (int m = 0; m < segGammai.length; m++) {
					
					double SUMMATION = 0.0;
					
					for (int j = 0; j < ncomps; j++) {
						if(pure && i!=j)
							continue; // if pure substance, skip all other
					
						Compound compj = comps[j];
						double []areaj = pure ? compj.area : PROFILE[j];
						double []segGammaj = segGamma[j];
						
						for(int n = 0; n < areaj.length; n++) {

							double deltaW_HB = 0;

							double sigma_mn = compi.sigmaAvg[m]+compj.sigmaAvg[n];
							double sigma_mnDelta = compi.sigmaDelta[m] + compj.sigmaDelta[n];
							sigma_mn = (sigma_mn + fcorr*sigma_mnDelta)/f_ortho;
							double deltaW = (fpol*alpha/2.0)* sigma_mn*sigma_mn;
							
							double expDeltaW_RT = Math.exp(-(deltaW + deltaW_HB) * inv_RT);
							
							SUMMATION += areaj[n] * segGammaj[n]*expDeltaW_RT;
						}
					}
					// substitute with the new value
					segGammai[m] = factor/SUMMATION;
				}
			}
			
			double newnorm = 0;
			for (int i = 0; i < ncomps; i++) {
				newnorm += blas.dnrm2(segGamma[i].length, segGamma[i], 1);
			}
			newnorm /= ncomps;
			
			if(Math.abs((norm - newnorm)/newnorm) <= tol || niter>maxiter)
				break;
			norm = newnorm;

			++niter;
//			System.out.println("SEGGAMMA niter:" + niter + " norm:" + norm);
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
			averageACOSMO += z[i]*comps[i].areaTotal;
		}

		for (int i = 0; i < ncomps; i++) {
			double []areai = comps[i].area;
			for(int m=0; m < areai.length; ++m){
				PROFILE[i][m] = z[i]*areai[m]/averageACOSMO;
			}
		}
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
		aEff = 7.5; // FIXME: remove this line
		alpha = 0.3/E0*Math.pow(aEff, 1.5);
		
		RNORM = new double[ncomps];
		QNORM = new double[ncomps];
		for(int i=0; i<ncomps; ++i){
			for (int m = 0; m < comps[i].area.length; m++) {
				// initialize all SEGGAMMAPR (reused between calculations)
				SEGGAMMAPR[i][m] = 1.0;
				SEGGAMMA[i][m] = 1.0;
			}
		}
		// some more composition independent properties
		for(int i=0; i<ncomps; ++i){
			RNORM[i] = comps[i].volumeTotal;
			QNORM[i] = comps[i].areaTotal/anorm;
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
		// Solve for the SEGGAMMA of the mixture
		solveGamma(false, SEGGAMMA, TOLERANCE);
		
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
		aEff = 7.5; // FIXME: remove this line

		for(int i=0; i<ncomps; ++i){
			double phi_zi = RNORM[i]/sum_zi_ri;
			double psi_zi = Math.pow(RNORM[i], rPower)/sum_ziR_ri;
			double theta_zi = QNORM[i]/sum_zi_qi;

			double lnGammaSG = 0;
			if(!ignoreSG){
				// expression from DOI:10.1021/ie070465z
				double phi_theta = phi_zi/theta_zi;
				lnGammaSG += Math.log(psi_zi) + 1 - psi_zi - (coord/2)*QNORM[i]*
					( Math.log(phi_theta) + 1 - phi_theta);
			}

			// CALCULATION OF GAMMAS
			double lnGammaRestoration = 0.0;
			double areai[] = comps[i].area;
			for(int m=0; m < areai.length; ++m){
				double lnMixSeg = Math.log(SEGGAMMA[i][m]);
				double lnMixPure = Math.log(SEGGAMMAPR[i][m]);
				
				lnGammaRestoration += (areai[m]/aEff)*(lnMixSeg - lnMixPure);
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
