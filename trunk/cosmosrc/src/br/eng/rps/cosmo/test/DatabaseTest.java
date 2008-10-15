package br.eng.rps.cosmo.test;

import junit.framework.TestCase;
import br.eng.rps.cosmo.COSMOSAC;
import br.eng.rps.cosmo.COSMOSACCompound;
import br.eng.rps.cosmo.COSMOSACDataBase;

public class DatabaseTest extends TestCase {

	public void testMethylAcetateWater() throws Exception{
		double T =330.15;

		COSMOSACDataBase db = COSMOSACDataBase.getInstance();
		COSMOSACCompound c1 = db.getComp("Methyl-Acetate");
		COSMOSACCompound c2 = db.getComp("Water");

		assertEquals(c1.charge.length, c1.sigma.length);
		assertEquals(c1.sigma.length, c2.sigma.length);

		double cavityVolume[] = new double[2];
		cavityVolume[0] = c1.Vcosmo;
		cavityVolume[1] = c2.Vcosmo;

		double [][] sigma = new double[2][];
		sigma[0] = c1.sigma;
		sigma[1] = c2.sigma;

		COSMOSAC cosmosac = new COSMOSAC();
		cosmosac.setParameters(cavityVolume, c1.charge, sigma);

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

	public void testChloroform14Dioxane() throws Exception{
		double T = 50 + 273.15;

		COSMOSACDataBase db = COSMOSACDataBase.getInstance();
		COSMOSACCompound c1 = db.getComp("chloroform");
		COSMOSACCompound c2 = db.getComp("1,4-dioxane");

		assertEquals(c1.charge.length, c1.sigma.length);
		assertEquals(c1.sigma.length, c2.sigma.length);

		double cavityVolume[] = new double[2];
		cavityVolume[0] = c1.Vcosmo;
		cavityVolume[1] = c2.Vcosmo;

		double [][] sigma = new double[2][];
		sigma[0] = c1.sigma;
		sigma[1] = c2.sigma;

		COSMOSAC cosmosac = new COSMOSAC();
		cosmosac.setParameters(cavityVolume, c1.charge, sigma);

		cosmosac.setTemperature(T);

		// testing several compositions
		double [] z = new double[2];
		double [] lnGamma = new double[2];

		z[0] = 0.0932; z[1] = 1-z[0];
		cosmosac.setComposition(z);
		cosmosac.activityCoefficient(lnGamma);
		assertEquals(Math.exp(-0.722), Math.exp(lnGamma[0]), 0.5);
		assertEquals(Math.exp(0.004), Math.exp(lnGamma[1]), 0.5);

		z[0] = 0.3615; z[1] = 1-z[0];
		cosmosac.setComposition(z);
		cosmosac.activityCoefficient(lnGamma);
		assertEquals(Math.exp(-0.486), Math.exp(lnGamma[0]), 0.5);
		assertEquals(Math.exp(-0.057), Math.exp(lnGamma[1]), 0.5);

		z[0] = 0.9398; z[1] = 1-z[0];
		cosmosac.setComposition(z);
		cosmosac.activityCoefficient(lnGamma);
		assertEquals(Math.exp(-0.002), Math.exp(lnGamma[0]), 0.5);
		assertEquals(Math.exp(-0.972), Math.exp(lnGamma[1]), 0.5);
	}


	/**
	 * Smith et al. Table 12.2
	 * @throws Exception
	 */
	public void testMethylEthylCetoneToluene() throws Exception{
		double T = 50 + 273.15;

		COSMOSACDataBase db = COSMOSACDataBase.getInstance();
		COSMOSACCompound c1 = db.getComp("Methyl-Ethyl-Ketone");
		COSMOSACCompound c2 = db.getComp("Toluene");

		assertEquals(c1.charge.length, c1.sigma.length);
		assertEquals(c1.sigma.length, c2.sigma.length);

		double cavityVolume[] = new double[2];
		cavityVolume[0] = c1.Vcosmo;
		cavityVolume[1] = c2.Vcosmo;

		double [][] sigma = new double[2][];
		sigma[0] = c1.sigma;
		sigma[1] = c2.sigma;

		COSMOSAC cosmosac = new COSMOSAC();
		cosmosac.setParameters(cavityVolume, c1.charge, sigma);

		cosmosac.setTemperature(T);

		// testing several compositions
		double [] z = new double[2];
		double [] lnGamma = new double[2];

		z[0] = 0.0895; z[1] = 1-z[0];
		cosmosac.setComposition(z);
		cosmosac.activityCoefficient(lnGamma);
		assertEquals(Math.exp(0.266), Math.exp(lnGamma[0]), 0.5);
		assertEquals(Math.exp(0.009), Math.exp(lnGamma[1]), 0.5);

		z[0] = 0.5119; z[1] = 1-z[0];
		cosmosac.setComposition(z);
		cosmosac.activityCoefficient(lnGamma);
		assertEquals(Math.exp(0.043), Math.exp(lnGamma[0]), 0.5);
		assertEquals(Math.exp(0.100), Math.exp(lnGamma[1]), 0.5);

		z[0] = 0.9102; z[1] = 1-z[0];
		cosmosac.setComposition(z);
		cosmosac.activityCoefficient(lnGamma);
		assertEquals(Math.exp(-0.003), Math.exp(lnGamma[0]), 0.5);
		assertEquals(Math.exp(0.237), Math.exp(lnGamma[1]), 0.5);
	}

	public void testBenzene224TrimethylPentane() throws Exception{
		double T = 50 + 273.15;

		COSMOSACDataBase db = COSMOSACDataBase.getInstance();
		COSMOSACCompound c1 = db.getComp("Benzene");
		COSMOSACCompound c2 = db.getComp("2,2,4-trimethylpentane");

		assertEquals(c1.charge.length, c1.sigma.length);
		assertEquals(c1.sigma.length, c2.sigma.length);

		double cavityVolume[] = new double[2];
		cavityVolume[0] = c1.Vcosmo;
		cavityVolume[1] = c2.Vcosmo;

		double [][] sigma = new double[2][];
		sigma[0] = c1.sigma;
		sigma[1] = c2.sigma;

		COSMOSAC cosmosac = new COSMOSAC();
		cosmosac.setParameters(cavityVolume, c1.charge, sigma);

		cosmosac.setTemperature(T);

		// testing several compositions
		double [] z = new double[2];
		double [] lnGamma = new double[2];

		z[0] = 0.4; z[1] = 1-z[0];
		cosmosac.setComposition(z);
		cosmosac.activityCoefficient(lnGamma);
		assertEquals(Math.exp(0.23), Math.exp(lnGamma[0]), 0.5);
		assertEquals(Math.exp(0.07), Math.exp(lnGamma[1]), 0.5);

		z[0] = 0.6; z[1] = 1-z[0];
		cosmosac.setComposition(z);
		cosmosac.activityCoefficient(lnGamma);
		assertEquals(Math.exp(0.16), Math.exp(lnGamma[0]), 0.5);
		assertEquals(Math.exp(0.18), Math.exp(lnGamma[1]), 0.5);

		z[0] = 0.9; z[1] = 1-z[0];
		cosmosac.setComposition(z);
		cosmosac.activityCoefficient(lnGamma);
		assertEquals(Math.exp(0.01), Math.exp(lnGamma[0]), 0.5);
		assertEquals(Math.exp(0.49), Math.exp(lnGamma[1]), 0.5);
	}
}