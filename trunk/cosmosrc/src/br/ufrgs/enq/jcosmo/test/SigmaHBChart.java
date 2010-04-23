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
public class SigmaHBChart {

	public static void main(String[] args) throws Exception {
		
		JFrame dlg = new JFrame();
		dlg.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		dlg.setLayout(new BorderLayout());
		
		String name = "ACETONE";
//		String name = "METHANE";
//		String name = "ETHANOL";
//		String name = "METHANOL";
//		String name = "WATER";

		String folder = "moltest/";
		String extension = ".pcm.gout";
		
		SigmaProfilePanel chart = new SigmaProfilePanel(name);
		dlg.add(chart, BorderLayout.CENTER);
		
		String fileName = folder + name + extension;		
		SigmaProfileGenerator sigmaParser;
		sigmaParser = new SigmaProfileGenerator(SigmaProfileGenerator.FileType.GAMESS_PCM);
		
		
		sigmaParser.parseFile(fileName, SigmaProfileGenerator.RAV);
		
		double sigma[] = sigmaParser.getAveragedChargeDensity();
		double area[] = sigmaParser.getOriginalArea();
		int[] elem = sigmaParser.getElem();
		
		sigmaParser.simpleSorting(area, sigma);
		chart.addProfile("Sigma", sigmaParser.getChargeDensity(), sigmaParser.getSortedArea());
		
		double[] areaHB = new double[area.length];
		
		for (int m = 0; m < area.length; m++) {
			if(elem[m]==1 || elem[m]==7 || elem[m]==8){
				areaHB[m] = area[m];
				area[m] = 0;
			}
		}

		sigmaParser.simpleSorting(area, sigma);
		chart.addProfile("Sigma No HB", sigmaParser.getChargeDensity(), sigmaParser.getSortedArea());
		
		sigmaParser.simpleSorting(areaHB, sigma);
		chart.addProfile("Sigma HB", sigmaParser.getChargeDensity(), sigmaParser.getSortedArea());
		
		dlg.setSize(600, 400);
		dlg.setLocationRelativeTo(null);
		dlg.setVisible(true);
		
	}
}
