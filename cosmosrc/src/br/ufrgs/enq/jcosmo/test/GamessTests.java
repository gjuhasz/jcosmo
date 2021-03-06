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

import br.ufrgs.enq.jcosmo.COSMOSAC;
import br.ufrgs.enq.jcosmo.COSMOSACCompound;
import br.ufrgs.enq.jcosmo.COSMOSACDataBase;
import br.ufrgs.enq.jcosmo.COSMOSAC_G;

public class GamessTests {

	public static void main(String[] args) throws Exception {
		
		COSMOSACDataBase db = COSMOSACDataBase.getInstance();
		
		COSMOSACCompound comps[] = new COSMOSACCompound[2];
		comps[0] = db.getComp("ethyl-acetate");
		comps[1] = db.getComp("water");
		
		COSMOSAC cosmosac = new COSMOSAC_G();
//		COSMOSAC cosmosac = new COSMOSAC();
		
//		cosmosac.setSigmaHB(COSMOSAC.SIGMAHB);
		
		
		cosmosac.setComponents(comps);
		cosmosac.setTemperature(373.15);

		double x[] = new double[2];
		double lnGamma[] = new double[2];

		x[0] = 0.5;
		x[1] = 1-x[0];
			
		cosmosac.setComposition(x);
		cosmosac.activityCoefficient(lnGamma);
			
		System.out.println(x[0] +"\t"+ Math.exp(lnGamma[0]) +"\t"+ Math.exp(lnGamma[1]) +"\t"+ lnGamma[0] +"\t"+ lnGamma[1]);
	}
}
