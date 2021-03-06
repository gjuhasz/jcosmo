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

import br.ufrgs.enq.jcosmo.COSMOPACMulti;
import br.ufrgs.enq.jcosmo.COSMOSACCompound;

/**
 * This program mimics the original FORTRAN code from VT-2005
 * @author rafael
 *
 */
public class COSMOPACMultiTest {

	public static void main(String[] args) {

		double T =330.15;
		
		COSMOPACMulti cosmosac = new COSMOPACMulti();
		
		COSMOSACCompound comps[] = new COSMOSACCompound[2];
		comps[0] = new COSMOSACCompound();
		comps[0].name = "METHYL-ACETATE";
		comps[1] = new COSMOSACCompound();
		comps[1].name = "WATER";
			
		try {
			cosmosac.setComponents(comps);
		} catch (Exception e) {
			e.printStackTrace();
		}
		
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
