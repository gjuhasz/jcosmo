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
import br.ufrgs.enq.jcosmo.SigmaProfileGenerator;

/**
 * Test for multiple sigma descriptors.
 * 
 * @author rafael
 *
 */
public class SigmaAveragingTest {

	public static void main(String[] args) throws Exception {
		
//		String name = "ACETONE";
//		String name = "METHANE";
		String name = "ETHANOL";
//		String name = "METHANOL";
//		String name = "DIMETHYL_SULFOXIDE";
		
		String folder = "profiles/RM1/";
		String extension = ".cos";
		
		String fileName = folder + name + extension;
		
		SigmaProfileGenerator sigmaParser;
		
		sigmaParser = new SigmaProfileGenerator(SigmaProfileGenerator.FileType.MOPAC);
		sigmaParser.parseFile(fileName, COSMOSAC.RAV);
		
		double[] area = sigmaParser.getOriginalArea();
		double []sigma = sigmaParser.getOriginalChargeDensity();
		double[] sigmaAvg = sigmaParser.getAveragedChargeDensity();
		
		double sigmaT = 0, sigmaAvgT = 0;
		for (int i = 0; i < sigmaAvg.length; i++) {
			sigmaT += area[i]*sigma[i];
			sigmaAvgT += area[i]*sigmaAvg[i];
		}
		
		System.out.println("Total charge: " + sigmaT);
		System.out.println("Average total charge: " + sigmaAvgT);
	}
}
