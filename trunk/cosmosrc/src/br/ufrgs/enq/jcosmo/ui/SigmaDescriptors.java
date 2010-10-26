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
import java.io.File;
import java.text.DecimalFormat;
import java.text.NumberFormat;

import javax.swing.JButton;
import javax.swing.JComboBox;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JTextField;

import br.ufrgs.enq.jcosmo.COSMOPAC;
import br.ufrgs.enq.jcosmo.COSMOSAC;
import br.ufrgs.enq.jcosmo.COSMOSAC_G;
import br.ufrgs.enq.jcosmo.MolParser;
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

		final JTextField nameField = new JTextField("HYDROGEN_FLUORIDE", 16);
		nameField.setText("N-FORMYLMORPHOLINE");
		String []fileTypeList = {"MOPAC", "SVP-GAMESS", "PCM-GAMESS", "COSMO-GAMESS"};
		final JComboBox fileType = new JComboBox(fileTypeList);
		String []analysisTypeList = {"Atom Type", "None", "Polarizability", "Low", "H, N, O, ..."};
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
				
				NumberFormat nf = NumberFormat.getInstance();
				nf.setMaximumFractionDigits(2);

				String folder = "moltest/";
				FileType type = SigmaProfileGenerator.FileType.MOPAC;
				String extension = ".cos";
				COSMOSAC model = null;
				double sigmaHBmin = 0.0038;

				if(fileType.getSelectedItem().equals("MOPAC")){
					model = new COSMOPAC();
//					folder = "mopac/";
//					folder = "mopPOA1_all/";
//					folder = "mopRM1/";
//					folder = "profiles/POA1/";
					folder = "profiles/RM1_1.18/";
					folder = "profiles/RM1/";
				}
				else if(fileType.getSelectedItem().equals("SVP-GAMESS")){
					folder = "moltest/";
					extension = ".svp.gout";
					type = SigmaProfileGenerator.FileType.GAMESS_SVP;
					model = new COSMOSAC();
				}
				else if(fileType.getSelectedItem().equals("PCM-GAMESS")){
					folder = "mopac/";
					extension = ".pcm.gout";
					type = SigmaProfileGenerator.FileType.GAMESS_PCM;
					model = new PCMSAC();
				}
				else if(fileType.getSelectedItem().equals("COSMO-GAMESS")){
//					folder = "moltest/";
//					folder = "gam6-31++G2d,p/";
//					folder = "gamSTO3/";
//					folder = "gam6-311G/";
//					folder = "gam6-31Gd/";
					folder = "profiles/gamess/";
					extension = ".gout";
					type = SigmaProfileGenerator.FileType.GAMESS;
					model = new COSMOSAC_G();
				}

				double rav = model.getRav();

				String fileName = folder + nameField.getText() + extension;

				SigmaProfileGenerator sigmaParser = new SigmaProfileGenerator(type);

				try{
					if(analysisType.getSelectedItem().equals("None")){
						sigmaParser.parseFile(fileName, rav);
						
						double[] sigmaBase = sigmaParser.getAveragedChargeDensity();
						double[] area = sigmaParser.getOriginalArea();
						double areaT = 0;
						double sigmaAvg = 0;
						for (int i = 0; i < area.length; i++) {
							areaT += area[i];
							sigmaAvg += area[i]*Math.abs(sigmaBase[i]);
						}
						sigmaAvg /= areaT;

					
						chart.addProfile("Sigma (" + areaT + ") Abs. Avg*1000=" + sigmaAvg*1000, sigmaParser.getChargeDensity(), sigmaParser.getSortedArea());
					}

					if(analysisType.getSelectedItem().equals("Atom Type")){
						sigmaParser.parseFile(fileName, rav);
						double[] sigmaBase = sigmaParser.getAveragedChargeDensity();

						double[] area = sigmaParser.getOriginalArea();
						int[] atoms = sigmaParser.getAtom();
						
						MolParser molParser = new MolParser();
						String name = nameField.getText();
						File molFile = new File(folder + name + ".mol");
						molParser.parseFile(molFile.getPath());
						
						// lets filter the H-[N,O,F,Cl,Br,I] atoms, the donnors
						double[] areaHDonnor = new double[area.length];
						for (int m = 0; m < area.length; m++) {
							int boundedType[] = {7, 8, 9, 17, 35, 53};
							if(sigmaBase[m]<-sigmaHBmin && molParser.matchType(atoms[m], 1, boundedType)){
								areaHDonnor[m] += area[m];
								area[m] = 0;
							}
						}

						// lets filter the [N,O,...]-H atoms (HB-acceptor bounded to H)
						double[] areaNH = new double[area.length];
						for (int m = 0; m < area.length; m++) {
							int atomType[] = {7, 8, 9, 17, 35, 53};
							if(sigmaBase[m]>sigmaHBmin && molParser.matchType(atoms[m], atomType, 1)){
								areaNH[m] += area[m];
								area[m] = 0;
							}
						}
						
						// lets filter the [N, O, ...] atoms
						double[] areaN = new double[area.length];
						for (int m = 0; m < area.length; m++) {
							int atomType[] = {7, 8, 9, 17, 35, 53};
//							All on group2
//							if(sigma1[m]>0 && molParser.matchType(atoms[m], atomType, 0)){
//								area2[m] += area0[m];
//								area0[m] = 0;
//							}

							// Oxygen goes to area2 only if double bounded
							if(sigmaBase[m]>sigmaHBmin && molParser.matchBondType(atoms[m], atomType, 2)){
								areaN[m] += area[m];
								area[m] = 0;
							}
							if(sigmaBase[m]>sigmaHBmin && molParser.matchBondType(atoms[m], atomType, 3)){
								areaN[m] += area[m];
								area[m] = 0;
							}
							// single bounded Oxygen go to area1
							// FIXME: detect automatically ACETATE and OXIDE oxygens
							if(!nameField.getText().endsWith("ATE") && !nameField.getText().endsWith("OXIDE")){
								if(sigmaBase[m]>sigmaHBmin && molParser.matchType(atoms[m], 8, 0)){
									areaNH[m] += area[m];
									area[m] = 0;
								}
							}
							if(sigmaBase[m]>sigmaHBmin && molParser.matchType(atoms[m], 8, 0)){
								areaNH[m] += area[m];
								area[m] = 0;
							}
							if(sigmaBase[m]>sigmaHBmin && molParser.matchType(atoms[m], atomType, 0)){
								areaN[m] += area[m];
								area[m] = 0;
							}
						}
						
						// lets filter the H-C-[N,O,F,Cl,Br,I] atoms, the donnors 2
						double[] areaHDonnor2 = new double[area.length];
						for (int m = 0; m < area.length; m++) {
							int boundedType2[] = {7, 8, 9, 17, 35, 53};
							if(sigmaBase[m]<sigmaHBmin && molParser.matchType(atoms[m], 1, 6, boundedType2)){
								areaHDonnor2[m] += area[m];
								area[m] = 0;
							}
						}

						
						double partial;

						// adding all other elements
						partial = sigmaParser.simpleSorting(area, sigmaBase);
						chart.addProfile("0. NonHB (" + nf.format(partial) + ")", sigmaParser.getChargeDensity(), sigmaParser.getSortedArea());
						
						partial = sigmaParser.simpleSorting(areaHDonnor, sigmaBase);
						chart.addProfile("1. H-[N, O, ...] (" + nf.format(partial) + ")", sigmaParser.getChargeDensity(), sigmaParser.getSortedArea());
						
						partial = sigmaParser.simpleSorting(areaNH, sigmaBase);
						chart.addProfile("2. [N, O, ...]-H (" + nf.format(partial) + ")", sigmaParser.getChargeDensity(), sigmaParser.getSortedArea());
						
						partial = sigmaParser.simpleSorting(areaN, sigmaBase);
						chart.addProfile("3. [N, O, ...] (" + nf.format(partial) + ")", sigmaParser.getChargeDensity(), sigmaParser.getSortedArea());
						
						partial = sigmaParser.simpleSorting(areaHDonnor2, sigmaBase);
						chart.addProfile("4. H-C-[N, O, ...] (" + nf.format(partial) + ")", sigmaParser.getChargeDensity(), sigmaParser.getSortedArea());
					}

					// Polarizability analysis
					if(analysisType.getSelectedItem().equals("Polarizability")){
						
						rav = 1.0;
						double fortho = 0.79209;
						double rav2 = 1.5*rav;
						
						sigmaParser.parseFile(fileName, rav);
						double[] sigmaBase = sigmaParser.getAveragedChargeDensity();

						sigmaParser.parseFile(fileName, rav2);
						double[] sigma2 = sigmaParser.getAveragedChargeDensity();

						double[] area = sigmaParser.getOriginalArea();
						double[] sigmaT = new double[area.length];

						sigmaParser.simpleSorting(area, sigmaBase);

						for (int m = 0; m < area.length; m++) {
//							sigmaT[m] = 1000*Math.abs(fcorr*sigmaBase[m] - sigma2[m]);
							sigmaT[m] = 1000*(sigma2[m] - fortho*sigmaBase[m]);
						}
//						double []sT = {-2, 0, 2};
						double []sT = {-1, -0.5, 0.5, 1};
//						double []sT = {0.7, 2, 3};
//						double []sT = {1, 2, 2.5};
//						double []sT = {1.5, 2};
//						double []sT = {-1, 0.18, 0.18, 1};
						double[] areaT = new double[area.length];
						for (int i = -1; i < sT.length; i++) {
							double stLow, stUp;
							double partialArea = 0;
							if(i<0){
								stLow = -Double.MAX_VALUE;
								stUp = sT[0];
							}
							else if(i<sT.length-1){
								stLow = sT[i];
								stUp = sT[i+1];
							}
							else{
								stLow = sT[i];
								stUp = Double.MAX_VALUE;
							}

							for (int m = 0; m < area.length; m++) {
								if(sigmaT[m]>=stLow && sigmaT[m]<stUp){
									areaT[m] = area[m];
									partialArea += area[m];
								}
								else
									areaT[m] = 0;
							}
							sigmaParser.simpleSorting(areaT, sigmaBase);
							
							DecimalFormat fm = new DecimalFormat();
							fm.setMaximumFractionDigits(2);
							String partial = String.format(" (%s)", fm.format(partialArea));
							if(i==-1)
								chart.addProfile("sT<" + stUp + partial, sigmaParser.getChargeDensity(), sigmaParser.getSortedArea());
							else if(i<sT.length-1)
								chart.addProfile("sT=" + stLow + " to " + stUp + partial, sigmaParser.getChargeDensity(), sigmaParser.getSortedArea());
							else
								chart.addProfile("sT>" + stLow + partial, sigmaParser.getChargeDensity(), sigmaParser.getSortedArea());
						}
					}
					if(analysisType.getSelectedItem().equals("Low")){
						sigmaParser.parseFile(fileName, rav);
						double[] sigmaBase = sigmaParser.getAveragedChargeDensity();

						fileName = folder + nameField.getText() + ".low" + extension;
						sigmaParser.parseFile(fileName, rav);
						double[] sigma2 = sigmaParser.getAveragedChargeDensity();

						// value from Klamt (COSMO-RS refinement)
						double fcorr = 0.816;
//						fcorr = 0.4;

						double[] area = sigmaParser.getOriginalArea();
						double[] sigmaT = new double[area.length];

						sigmaParser.simpleSorting(area, sigmaBase);

						for (int m = 0; m < area.length; m++) {
//							sigmaT[m] = 1000*(fcorr*sigma1[m] - sigma2[m]);
							sigmaT[m] = 1000*Math.abs(sigma2[m]-fcorr*sigmaBase[m]);
//							sigmaT[m] = 10*(Math.abs(sigmaBase[m])/(Math.abs(sigma2[m]) + 1e-5)-1.38);
//							sigmaT[m] = 10*Math.abs(sigma2[m]-fcorr*sigmaBase[m]);
						}
//						double []sT = {-2, 0, 2};
//						double []sT = {0.4, 1};
						double []sT = {0.5, 1};
//						double []sT = {0.7, 2, 3};
//						double []sT = {1, 2, 2.5};
//						double []sT = {1.5, 2};
//						double []sT = {-1, 0.18, 0.18, 1};
						double[] areaT = new double[area.length];
						for (int i = -1; i < sT.length; i++) {
							double stLow, stUp;
							double partialArea = 0;
							if(i<0){
								stLow = -Double.MAX_VALUE;
								stUp = sT[0];
							}
							else if(i<sT.length-1){
								stLow = sT[i];
								stUp = sT[i+1];
							}
							else{
								stLow = sT[i];
								stUp = Double.MAX_VALUE;
							}

							for (int m = 0; m < area.length; m++) {
								if(sigmaT[m]>=stLow && sigmaT[m]<stUp){
									areaT[m] = area[m];
									partialArea += area[m];
								}
								else
									areaT[m] = 0;
							}
							sigmaParser.simpleSorting(areaT, sigmaBase);
							
							DecimalFormat fm = new DecimalFormat();
							fm.setMaximumFractionDigits(2);
							String partial = String.format(" (%s)", fm.format(partialArea));
							if(i==-1)
								chart.addProfile("sT<" + stUp + partial, sigmaParser.getChargeDensity(), sigmaParser.getSortedArea());
							else if(i<sT.length-1)
								chart.addProfile("sT=" + stLow + " to " + stUp + partial, sigmaParser.getChargeDensity(), sigmaParser.getSortedArea());
							else
								chart.addProfile("sT>" + stLow + partial, sigmaParser.getChargeDensity(), sigmaParser.getSortedArea());
						}
					}
					if(analysisType.getSelectedItem().equals("H, N, O, ...")){
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
