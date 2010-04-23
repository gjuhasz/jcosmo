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

import br.ufrgs.enq.jcosmo.SigmaProfileGenerator;

/**
 * Test for multiple sigma descriptors.
 * 
 * @author rafael
 *
 */
public class SigmaDescriptors {

	public static void main(String[] args) throws Exception {
		
		String name = "WATER";
//		String name = "METHANE";
		
		SigmaProfileGenerator sigmaParser;
		
		sigmaParser = new SigmaProfileGenerator(SigmaProfileGenerator.FileType.GAMESS,
				SigmaProfileGenerator.RAV, 51);
		sigmaParser.parseFile("moltest/" + name + ".gout");
		
		double[] area = sigmaParser.getOriginalArea();
		
		double[] sigma1 = sigmaParser.getAveragedChargeDensity();
		
		sigmaParser = new SigmaProfileGenerator(SigmaProfileGenerator.FileType.GAMESS,
				SigmaProfileGenerator.RAV*2, 51);
		sigmaParser.parseFile("moltest/" + name + ".gout");
		
		double[] sigma2 = sigmaParser.getAveragedChargeDensity();
		
		double fcorr = 0.816;
//		double fcorr = COSMOSAC.FPOL;
		
		System.out.println("Sigma(RAV), Sigma(RAV*2), SigmaT");
		for (int i = 0; i < area.length; i++) {
			System.out.println(sigma1[i] + ", " + sigma2[i] + ", " + 1000*(sigma2[i]-fcorr*sigma1[i]));
		}
	}
}
