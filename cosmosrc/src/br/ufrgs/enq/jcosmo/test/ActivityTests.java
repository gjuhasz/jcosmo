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

import java.sql.SQLException;

import br.ufrgs.enq.jcosmo.COSMOSAC;
import br.ufrgs.enq.jcosmo.COSMOSACCompound;
import br.ufrgs.enq.jcosmo.COSMOSACDataBase;

public class ActivityTests {

	public static void main(String[] args) {
		
		COSMOSACDataBase db = COSMOSACDataBase.getInstance();
		COSMOSACCompound c1;
		COSMOSACCompound c2;
		try {
			c1 = db.getComp("methyl-acetate");
			c2 = db.getComp("water");
		} catch (SQLException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			return;
		}

		double[] cavityVolume = new double[2];
		cavityVolume[0] = c1.Vcosmo;
		cavityVolume[1] = c2.Vcosmo;

		double [][] sigma = new double[2][];
		sigma[0] = c1.sigma;
		sigma[1] = c2.sigma;

		COSMOSAC cosmosac = new COSMOSAC();
		cosmosac.setParameters(cavityVolume, c1.charge, sigma);

		cosmosac.setTemperature(330.15);

		int n = 100;
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
