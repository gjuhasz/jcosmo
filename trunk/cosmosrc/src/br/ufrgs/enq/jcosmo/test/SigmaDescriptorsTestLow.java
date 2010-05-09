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
import br.ufrgs.enq.jcosmo.SigmaProfileGenerator.FileType;

/**
 * Test for multiple sigma descriptors.
 * 
 * @author rafael
 *
 */
public class SigmaDescriptorsTestLow {

	public static void main(String[] args) throws Exception {
		
//		String name = "ACETONE";
		String name = "METHANE";
//		String name = "ETHANOL";
//		String name = "METHANOL";
//		String name = "WATER";
//		String name = "CHLOROFORM";
//		String name = "CARBON_TETRACHLORIDE";
		
		String folder = "mopac/";
//		String extension = ".pcm.gout";
//		FileType type = SigmaProfileGenerator.FileType.GAMESS_PCM;
		String extension = ".gout";
		FileType type = SigmaProfileGenerator.FileType.GAMESS;
		
		String fileName;
		
		SigmaProfileGenerator sigmaParser;
		
		sigmaParser = new SigmaProfileGenerator(type);
		
		fileName = folder + name + extension;
		sigmaParser.parseFile(fileName, COSMOSAC.RAV);
		double[] sigma1 = sigmaParser.getAveragedChargeDensity();
		
		fileName = folder + name + ".low" + extension;
		sigmaParser.parseFile(fileName, COSMOSAC.RAV);
		double[] sigma2 = sigmaParser.getAveragedChargeDensity();
		
		
		double fcorr = 0.816;
//		double fcorr = 1;
//		double epsi = 0.5;
//		double fcorr = 1 - (epsi - 1)/(epsi + 0.5);
		
		double[] area = sigmaParser.getOriginalArea();
		int[] elem = sigmaParser.getElem();
		int[] atom = sigmaParser.getAtom();

		System.out.println("Atom, Elemnt, Area, Sigma, SigmaLow, SigmaT");
		for (int i = 0; i < area.length; i++) {
			System.out.println(atom[i] + ", " + elem[i] + ", " + area[i] + ", " +
					sigma1[i] + ", " + sigma2[i] + ", " + 1000*(sigma2[i]-fcorr*sigma1[i]));
		}
	}
}
