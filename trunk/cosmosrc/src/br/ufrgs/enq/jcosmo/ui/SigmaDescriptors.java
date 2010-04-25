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

package br.ufrgs.enq.jcosmo.ui;

import java.awt.BorderLayout;
import java.awt.FlowLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.JButton;
import javax.swing.JComboBox;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JTextField;

import br.ufrgs.enq.jcosmo.COSMOPAC;
import br.ufrgs.enq.jcosmo.COSMOSAC;
import br.ufrgs.enq.jcosmo.COSMOSAC_G;
import br.ufrgs.enq.jcosmo.PCMSAC;
import br.ufrgs.enq.jcosmo.SigmaProfileGenerator;
import br.ufrgs.enq.jcosmo.SigmaProfileGenerator.FileType;

/**
 * Test for multiple sigma descriptors.
 * 
 * @author rafael
 *
 */
public class SigmaDescriptors {

	public static void main(String[] args) throws Exception {

		JFrame dlg = new JFrame();
		dlg.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		dlg.setLayout(new BorderLayout());

		final SigmaProfileAreaPanel chart = new SigmaProfileAreaPanel();
		dlg.add(chart, BorderLayout.CENTER);

		final JTextField nameField = new JTextField("METHANOL", 20);
		String []fileTypeList = {"MOPAC", "PCM-GAMESS", "COSMO-GAMESS"};
		final JComboBox fileType = new JComboBox(fileTypeList);
		String []analysisTypeList = {"Polarizability", "Hydrogen-Bond"};
		final JComboBox analysisType = new JComboBox(analysisTypeList);

		JButton run = new JButton("Refresh");

		JPanel top = new JPanel(new FlowLayout());
		dlg.add(top, BorderLayout.PAGE_START);
		top.add(new JLabel("Name:"));
		top.add(nameField);
		top.add(new JLabel("File Type:"));
		top.add(fileType);
		top.add(new JLabel("Analysis:"));
		top.add(analysisType);
		top.add(run);
		
		ActionListener action = new ActionListener(){
			public void actionPerformed(ActionEvent evt) {				
				chart.clearProfiles();

				nameField.selectAll();
				nameField.grabFocus();

				String folder = "mopac/";
				FileType type = SigmaProfileGenerator.FileType.MOPAC;
				String extension = ".cos";
				COSMOSAC model = null;

				if(fileType.getSelectedItem().equals("MOPAC")){
					model = new COSMOPAC();
				}
				else if(fileType.getSelectedItem().equals("PCM-GAMESS")){
					folder = "moltest/";
					extension = ".pcm.gout";
					type = SigmaProfileGenerator.FileType.GAMESS_PCM;
					model = new PCMSAC();
				}
				else if(fileType.getSelectedItem().equals("COSMO-GAMESS")){
					folder = "moltest/";
					extension = ".gout";
					type = SigmaProfileGenerator.FileType.GAMESS;
					model = new COSMOSAC_G();
				}

				double rav = model.getRav();

				String fileName = folder + nameField.getText() + extension;

				SigmaProfileGenerator sigmaParser = new SigmaProfileGenerator(type);

				try{
					// Polarizability analysis
					if(analysisType.getSelectedIndex()==0){
						sigmaParser.parseFile(fileName, rav);
						double[] sigma1 = sigmaParser.getAveragedChargeDensity();

						sigmaParser.parseFile(fileName, rav*2);
						double[] sigma2 = sigmaParser.getAveragedChargeDensity();

						// value from Klamt (COSMO-RS refinement)
						double fcorr = 0.816;

						double[] area = sigmaParser.getOriginalArea();
						double[] sigmaT = new double[area.length];

						sigmaParser.simpleSorting(area, sigma1);

						for (int m = 0; m < area.length; m++) {
							//					sigmaT[m] = 1000*(sigma2[m]-fcorr*sigma1[m]);
							sigmaT[m] = 1000*Math.abs(sigma2[m]-fcorr*sigma1[m]);
						}

						double []sT = {0, 1.5, 3};
						double[] areaT = new double[area.length];
						for (int i = 0; i < sT.length; i++) {
							double stLow = sT[i];
							double stUp = Double.MAX_VALUE;
							if(i<sT.length-1)
								stUp = sT[i+1];

							for (int m = 0; m < area.length; m++) {
								if(sigmaT[m]>=stLow && sigmaT[m]<stUp)
									areaT[m] = area[m];
								else
									areaT[m] = 0;
							}
							sigmaParser.simpleSorting(areaT, sigma1);
							if(i<sT.length-1)
								chart.addProfile("sT=" + stLow + " to " + stUp, sigmaParser.getChargeDensity(), sigmaParser.getSortedArea());
							else
								chart.addProfile("sT>" + stLow, sigmaParser.getChargeDensity(), sigmaParser.getSortedArea());
						}
					}
					else{
						// Hydrogen-Bond analysis
						sigmaParser.parseFile(fileName, rav);
						
						double sigma[] = sigmaParser.getAveragedChargeDensity();
						double area[] = sigmaParser.getOriginalArea();
						int[] elem = sigmaParser.getElem();
						
						double[] areaHB = new double[area.length];
						for (int m = 0; m < area.length; m++) {
							if(elem[m]==1 || elem[m]==7 || elem[m]==8 || elem[m]==17){
								areaHB[m] = area[m];
								area[m] = 0;
							}
						}
						
						sigmaParser.simpleSorting(area, sigma);
						chart.addProfile("Sigma Others", sigmaParser.getChargeDensity(), sigmaParser.getSortedArea());
						
						sigmaParser.simpleSorting(areaHB, sigma);
						chart.addProfile("Sigma H, N, O, Cl", sigmaParser.getChargeDensity(), sigmaParser.getSortedArea());
					}
				}
				catch (Exception e) {
					e.printStackTrace();
				}
			}
		};

		nameField.addActionListener(action);
		run.addActionListener(action);

		dlg.setSize(800, 600);
		dlg.setLocationRelativeTo(null);
		dlg.setVisible(true);
		nameField.selectAll();
		nameField.grabFocus();
	}
}
