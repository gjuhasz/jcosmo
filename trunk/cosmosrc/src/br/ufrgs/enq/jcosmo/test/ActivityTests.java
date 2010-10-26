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

import br.ufrgs.enq.jcosmo.COSMOPAC;
import br.ufrgs.enq.jcosmo.COSMOSAC;
import br.ufrgs.enq.jcosmo.COSMOSACCompound;

public class ActivityTests {

	public static void main(String[] args) throws Exception {
		
		COSMOSACCompound comps[] = new COSMOSACCompound[2];
		comps[0] = new COSMOSACCompound();
		comps[1] = new COSMOSACCompound();
		comps[0].name = "N,N-DIMETHYLFORMAMIDE";
		comps[1].name = "N-HEXANE";
		
		comps[0].name = "N-OCTANE";
		comps[1].name = "ETHYL ACETATE";
		
		comps[0].name = "N-PENTANE";
		comps[1].name = "METHYL ETHYL KETONE";
		
		COSMOSAC cosmosac = new COSMOPAC();
		cosmosac.setComponents(comps);

		cosmosac.setTemperature(333.15);

		int n = 21;
		double x[] = new double[2];
		double lnGamma[] = new double[2];

		System.out.println("X1\tGamma1\tGamma2\tlnGamma1\tlnGamma2");

		for (int i = 0; i < n; i++) {
			x[0] = (double)i/(n-1);
			x[1] = 1-x[0];
			
			cosmosac.setComposition(x);
			cosmosac.activityCoefficient(lnGamma);
			
			System.out.println(x[0] +"\t"+ Math.exp(lnGamma[0]) +"\t"+ Math.exp(lnGamma[1]) +"\t"+ lnGamma[0] +"\t"+ lnGamma[1]);
		}
	}
}
