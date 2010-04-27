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

package br.ufrgs.enq.jcosmo.test;

import junit.framework.TestCase;
import br.ufrgs.enq.jcosmo.COSMOPACM;
import br.ufrgs.enq.jcosmo.COSMOSACCompound;
import br.ufrgs.enq.jcosmo.COSMOSACDataBase;

public class DatabaseTest extends TestCase {

	public void testMethylAcetateWater() throws Exception{
		double T =330.15;

		COSMOSACDataBase db = COSMOSACDataBase.getInstance();
		COSMOSACCompound c[] = new COSMOSACCompound[2];
		c[0] = db.getComp("Methyl-Acetate");
		c[1] = db.getComp("Water");

		COSMOPACM cosmosac = new COSMOPACM();
		cosmosac.setComponents(c);

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
		COSMOSACCompound c[] = new COSMOSACCompound[2];
		c[0] = db.getComp("chloroform");
		c[1] = db.getComp("1,4-dioxane");

		COSMOPACM cosmosac = new COSMOPACM();
		cosmosac.setComponents(c);

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
		COSMOSACCompound c[] = new COSMOSACCompound[2];
		c[0] = db.getComp("Methyl-Ethyl-Ketone");
		c[1] = db.getComp("Toluene");

		COSMOPACM cosmosac = new COSMOPACM();
		cosmosac.setComponents(c);

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
		COSMOSACCompound c[] = new COSMOSACCompound[2];
		c[0] = db.getComp("Benzene");
		c[1] = db.getComp("2,2,4-trimethylpentane");

		COSMOPACM cosmosac = new COSMOPACM();
		cosmosac.setComponents(c);

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