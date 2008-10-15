package br.eng.rps.cosmo.test;

import java.util.ArrayList;

import junit.framework.TestCase;
import br.eng.rps.cosmo.COSMOSAC;

public class COSMOSACTest extends TestCase {

	public void testMethylAcetateWater() throws Exception{
		double T =330.15;
		
		ArrayList<Double> chargeArgument = new ArrayList<Double>();
		ArrayList<Double> sigma1 = new ArrayList<Double>();	
		ArrayList<Double> sigma2 = new ArrayList<Double>();	
		
		// methyl acetate (638)
		CosmoImport.readSigmaProfile(638, chargeArgument, sigma1);		
		// water (1076)
		CosmoImport.readSigmaProfile(1076, chargeArgument, sigma2);	
		
		assertEquals(sigma1.size(), sigma2.size());
		
		double cavityVolume[] = {97.00036, 25.73454};
		
		double [] charge = new double[chargeArgument.size()];
		for (int i = 0; i < charge.length; i++) {
			charge[i] = chargeArgument.get(i);
		}
		double [][] sigma = new double[2][];
		sigma[0] = new double[charge.length];
		sigma[1] = new double[charge.length];
		for (int i = 0; i < charge.length; i++) {
			sigma[0][i] = sigma1.get(i);
			sigma[1][i] = sigma2.get(i);
		}
		
		COSMOSAC cosmosac = new COSMOSAC();
		cosmosac.setParameters(cavityVolume, charge, sigma);
		
		cosmosac.setTemperature(T);
		
		// testing several compositions
		double [] z = new double[2];
		double [] lnGamma = new double[2];
		z[0] = 0.00;
		
		System.out.println("SYSTEMP  " + T + " KELVIN");
		System.out.println("MOLE FRAC        GAMMA1         GAMMA2         LNGAMMA1       LNGAMMA2");

		while(z[0] < 1.0001){
			z[1] = 1-z[0];
			
			cosmosac.setComposition(z);
			
			cosmosac.activityCoefficient(lnGamma);
			
			System.out.println(" " + z[0] + "\t" + Math.exp(lnGamma[0]) + "\t" + Math.exp(lnGamma[1]) +
					"\t" + lnGamma[0] + "\t" + lnGamma[1] );
			
			z[0] += 0.005;
		}
	}
}
