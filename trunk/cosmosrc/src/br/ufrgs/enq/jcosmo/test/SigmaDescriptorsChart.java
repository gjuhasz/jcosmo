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

import java.awt.BorderLayout;

import javax.swing.JFrame;

import br.ufrgs.enq.jcosmo.SigmaProfileGenerator;
import br.ufrgs.enq.jcosmo.ui.SigmaProfilePanel;

/**
 * Test for multiple sigma descriptors.
 * 
 * @author rafael
 *
 */
public class SigmaDescriptorsChart {

	public static void main(String[] args) throws Exception {
		
		JFrame dlg = new JFrame();
		dlg.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		dlg.setLayout(new BorderLayout());
		
//		String name = "ACETONE";
//		String name = "METHANE";
//		String name = "ETHANOL";
		String name = "METHANOL";
//		String name = "WATER";
		
		String folder = "moltest/";
		String extension = ".pcm.gout";
		
		SigmaProfilePanel chart = new SigmaProfilePanel(name);
		dlg.add(chart, BorderLayout.CENTER);
		
		String fileName = folder + name + extension;
		
		SigmaProfileGenerator sigmaParser;
		
		sigmaParser = new SigmaProfileGenerator(SigmaProfileGenerator.FileType.GAMESS_PCM);
		
		sigmaParser.parseFile(fileName, SigmaProfileGenerator.RAV);
		double[] sigma1 = sigmaParser.getAveragedChargeDensity();
		
		sigmaParser.parseFile(fileName, SigmaProfileGenerator.RAV*2);
		double[] sigma2 = sigmaParser.getAveragedChargeDensity();
		
		double fcorr = 0.816;
//		double fcorr = COSMOSAC.FPOL;

		
		double[] area = sigmaParser.getOriginalArea();
		double[] sigmaT = new double[area.length];
		
		sigmaParser.simpleSorting(area, sigma1);
		chart.addProfile("full profile", sigmaParser.getChargeDensity(), sigmaParser.getSortedArea(), true);

		System.out.println("Atom, Elemnt, Area, Sigma(RAV), Sigma(RAV*2), SigmaT");
		for (int m = 0; m < area.length; m++) {
//			sigmaT[m] = 1000*(sigma2[m]-fcorr*sigma1[m]);
			sigmaT[m] = 1000*Math.abs(sigma2[m]-fcorr*sigma1[m]);
		}
		
//		double []sT = {-4, -1, 0, 1, 4};
		double []sT = {0, 1, 2, 4};
		double[] areaT = new double[area.length];
		for (int i = 0; i < sT.length-1; i++) {
			
			for (int m = 0; m < area.length; m++) {
				if(sigmaT[m]>=sT[i] && sigmaT[m]<sT[i+1])
					areaT[m] = area[m];
				else
					areaT[m] = 0;
			}
			sigmaParser.simpleSorting(areaT, sigma1);
			chart.addProfile("sT=" + sT[i] + " to " + sT[i+1], sigmaParser.getChargeDensity(), sigmaParser.getSortedArea());
		}
		
		dlg.setSize(600, 400);
		dlg.setLocationRelativeTo(null);
		dlg.setVisible(true);
	}
}
