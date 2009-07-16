package br.ufrgs.enq.jcosmo.test;

import junit.framework.TestCase;
import br.ufrgs.enq.jcosmo.COSMOPAC;
import br.ufrgs.enq.jcosmo.COSMOSACCompound;

/**
 * 
 * @author rafael
 *
 */
public class Cholesterol extends TestCase {
	
	/**
	 * Assert equals but with a relative tolerance
	 * @param value the correct value
	 * @param aprox the approximate value
	 * @param tol the tolerance
	 */
	public static void assertEquals(double value, double aprox, double tol){
		System.out.println("Expected:" + value + " found:" + aprox);
		TestCase.assertEquals(value, aprox, value * tol);
	}

	COSMOPAC geModel;
	
	double T = 308.2;
	static final double R = 8.314390E+00;
	
	COSMOSACCompound comps [] = new COSMOSACCompound[2];

	@Override
	protected void setUp() throws Exception {
		comps[0] = new COSMOSACCompound();
		comps[0].name = "WATER";		

		geModel = new COSMOPAC();
//		geModel.setIgnoreSG(true);
	}

	public void testCholesterolSol() throws Exception {
		comps[1] = new COSMOSACCompound();
		comps[1].name = "CHOLESTEROL";

		double z[] = new double[comps.length];
		z[0] = 1; z[1] = 1-z[0];
		
		double T = 308.2;
		
		double gammaInf[] = new double[comps.length];
		
		geModel.setComponents(comps);

		// find water solubility at 308.2 K
		geModel.setTemperature(T);
		geModel.setComposition(z);
		geModel.activityCoefficientLn(gammaInf, 0);
		
		for (int i = 0; i < gammaInf.length; i++) {
			gammaInf[i] = Math.exp(gammaInf[i]);
		}
		
		// data from "Thermal properties of cholesterol and estimation of impurity by
		// differential scanning calorimetry" Journal of Scientific Instruments 1969 Series 2 Volume 2
		double h_fus = 7.12e3 * 4.1868; // 7.12kcal/mol
		double Tm = 422.4;
		
		// xi^l gamma_i^l = zi gamma_i^s psi_i
		// considering that the solid is pure
		// zi gamma_i^s = 1
		// a approximation for psi_i is
		// psi_i = exp(h_fus/(R Tm) * (T - Tm)/T)
		
		double psi = Math.exp(h_fus/(R * Tm) * (T - Tm)/T);
		
		double sol = psi/gammaInf[1] * 1000 / 18;
		
		// Aqueous solubility
		double sol_exp = 3.7e-8; // M (mol/l)
		assertEquals(sol_exp, sol, 1.0);
	}
	
	public void testCampesterolSol() throws Exception {
		comps[1] = new COSMOSACCompound();
		comps[1].name = "CAMPESTEROL";

		double z[] = new double[comps.length];
		z[0] = 1; z[1] = 1-z[0];
		
		double T = 308.2;
		
		double gammaInf[] = new double[comps.length];
		
		geModel.setComponents(comps);

		// find water solubility at 308.2 K
		geModel.setTemperature(T);
		geModel.setComposition(z);
		geModel.activityCoefficientLn(gammaInf, 0);
		
		for (int i = 0; i < gammaInf.length; i++) {
			gammaInf[i] = Math.exp(gammaInf[i]);
		}
		
		// data from "Thermal properties of cholesterol and estimation of impurity by
		// differential scanning calorimetry" Journal of Scientific Instruments 1969 Series 2 Volume 2
		double h_fus = 7.12e3 * 4.1868; // 7.12kcal/mol
		double Tm = 139.93 + 273.15; // SST from DOI: 10.1021/jf0624289
		
		// xi^l gamma_i^l = zi gamma_i^s psi_i
		// considering that the solid is pure
		// zi gamma_i^s = 1
		// a approximation for psi_i is
		// psi_i = exp(h_fus/(R Tm) * (T - Tm)/T)
		
		double psi = Math.exp(h_fus/(R * Tm) * (T - Tm)/T);
		
		double sol = psi/gammaInf[1] * 1000 / 18;
		
		// Aqueous solubility
		double sol_exp = 1.6e-8; // M (mol/l)
		assertEquals(sol_exp, sol, 1.0);
	}
	
	public void testSitostanolSol() throws Exception {
		comps[1] = new COSMOSACCompound();
		comps[1].name = "SITOSTANOL";

		double z[] = new double[comps.length];
		z[0] = 1; z[1] = 1-z[0];
		
		double T = 308.2;
		
		double gammaInf[] = new double[comps.length];
		
		geModel.setComponents(comps);

		// find water solubility at 308.2 K
		geModel.setTemperature(T);
		geModel.setComposition(z);
		geModel.activityCoefficientLn(gammaInf, 0);
		
		for (int i = 0; i < gammaInf.length; i++) {
			gammaInf[i] = Math.exp(gammaInf[i]);
		}
		
		// data from "Thermal properties of cholesterol and estimation of impurity by
		// differential scanning calorimetry" Journal of Scientific Instruments 1969 Series 2 Volume 2
		double h_fus = 7.12e3 * 4.1868; // 7.12kcal/mol T = 422.4 K;
		double Tm = (273.15 + 144.76); // SSN from DOI: 10.1021/jf0624289
		
		// xi^l gamma_i^l = zi gamma_i^s psi_i
		// considering that the solid is pure
		// zi gamma_i^s = 1
		// a approximation for psi_i is
		// psi_i = exp(h_fus/(R Tm) * (T - Tm)/T)
		
		double psi = Math.exp(h_fus/(R * Tm) * (T - Tm)/T);
		
		double sol = psi/gammaInf[1] * 1000 / 18;
		
		// Aqueous solubility
		double sol_exp = 5.8e-9; // M (mol/l)
		assertEquals(sol_exp, sol, 1.0);
	}
	
	public void testSitosterolSol() throws Exception {
		comps[1] = new COSMOSACCompound();
		comps[1].name = "SITOSTEROL";

		double z[] = new double[comps.length];
		z[0] = 1; z[1] = 1-z[0];
		
		double T = 308.2;
		
		double gammaInf[] = new double[comps.length];
		
		geModel.setComponents(comps);

		// find water solubility at 308.2 K
		geModel.setTemperature(T);
		geModel.setComposition(z);
		geModel.activityCoefficientLn(gammaInf, 0);
		
		for (int i = 0; i < gammaInf.length; i++) {
			gammaInf[i] = Math.exp(gammaInf[i]);
		}
		
		// data from "Thermal properties of cholesterol and estimation of impurity by
		// differential scanning calorimetry" Journal of Scientific Instruments 1969 Series 2 Volume 2
		double h_fus = 7.12e3 * 4.1868; // 7.12kcal/mol
		double Tm = 138.5 + 273.15; // WSTSN from DOI: 10.1021/jf0624289
		
		// xi^l gamma_i^l = zi gamma_i^s psi_i
		// considering that the solid is pure
		// zi gamma_i^s = 1
		// a approximation for psi_i is
		// psi_i = exp(h_fus/(R Tm) * (T - Tm)/T)
		
		double psi = Math.exp(h_fus/(R * Tm) * (T - Tm)/T);
		
		double sol = psi/gammaInf[1] * 1000 / 18;
		
		// Aqueous solubility
		double sol_exp = 2.8e-9; // M (mol/l)
		assertEquals(sol_exp, sol, 1.0);
	}
	

	
}
