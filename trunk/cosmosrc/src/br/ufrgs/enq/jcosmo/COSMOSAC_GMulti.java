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

package br.ufrgs.enq.jcosmo;

import java.io.FileNotFoundException;
import java.util.HashMap;



/**
 * COSMO-SAC (GAMESS sigma profiles) activity model.
 *
 * @author Rafael de Pelegrini Soares
 * 
 */
public class COSMOSAC_GMulti extends COSMOSACMulti {
	private HashMap<String, COSMOSACCompound> compList;

	public String toString(){
		return "COSMO-SAC(GAMESS)";
	}
	
	public COSMOSAC_GMulti(){
		super(51, 3);
		
		// NONAQUEOUS, COST:0.2736488469786023 NP:167
//		idac/Alcohol-Alkane.csv AARD:0.3239770780689189 NP:10
//		idac/Alcohol-CycloAlkane.csv AARD:0.31659960715621416 NP:15
//		idac/Alkane-Alcohol.csv AARD:0.12249019700050272 NP:9
//		idac/Alkane-Alkane.csv AARD:0.021363199212819123 NP:3
//		idac/Alkane-AlkylHalide.csv AARD:0.28308602663247073 NP:5
//		idac/Alkane-Amine.csv AARD:0.31007765605690885 NP:4
//		idac/Alkane-CarboxilicAcid.csv AARD:NaN NP:0
//		idac/Alkane-Ketone.csv AARD:0.2026584124991165 NP:18
//		idac/Alkane-Phenol.csv AARD:0.13640569891573637 NP:22
//		idac/Alkene-Amine.csv AARD:0.22498699148647958 NP:3
//		idac/AlkylHalide-Alkane.csv AARD:0.3535428023306032 NP:13
//		idac/Amine-Alkane.csv AARD:0.6165306940215931 NP:5
//		idac/Aromatic-Alkane.csv AARD:0.1232432229575623 NP:6
//		idac/CycloAlkane-Alcohol.csv AARD:0.1432363420945877 NP:6
//		idac/CycloAlkane-AlkylHalide.csv AARD:0.029828597244639533 NP:5
//		idac/CycloAlkane-Amine.csv AARD:0.18564250391234083 NP:1
//		idac/CycloAlkane-CarboxilicAcid.csv AARD:0.7034224039183667 NP:2
//		idac/CycloAlkane-Phenol.csv AARD:0.14393537997518346 NP:5
//		idac/CarboxilicAcid-Alkane.csv AARD:2.558912369738963 NP:2
//		idac/CarboxilicAcid-CycloAlkane.csv AARD:2.128247009122318 NP:2
//		idac/CycloAlkene-Amine.csv AARD:0.010337077996353564 NP:3
//		idac/Ketone-Alcohol.csv AARD:0.2906249588540736 NP:14
//		idac/Ketone-Alkane.csv AARD:0.09654562496685164 NP:14
		setBeta(1.704710195);
		setFpol(0.49456570084204576);
		setFpol(1, 0.11917051421257);
		setFpol(2, 0.57751444534);
		setSigmaHB(0.0048793);
		setSigmaHB2(0.0);
		setSigmaHB3(1.0);
		setCHB(0);
		setCHB(2, 12318.276);
		setIgnoreSG(false);
		setCoord(10);
		setAnorm(107.16496305);
		setVnorm(66.69);
		
		// COST:0.44667417604299436 NP:206
		setBeta(0.8234057258502678);
		setBeta(1, 0.8234057258502678);
		setBeta(2, 0.8234057258502678);
		setCHB(0.0);
		setCHB(1, 0.0);
		setCHB(2, 11658.463520632726);
		setSigmaHB(0.0080);
		setSigmaHB2(0.0);
		setFpol(1.590397183396032);
		setFpol(1, 0.26171380795455323);
		setFpol(2, 0.43386856060593726);
		setIgnoreSG(false);
		setCoord(10.0);
		setAnorm(107.16496305);
		setVnorm(66.69);
	}
	
	public void setComponents(COSMOSACCompound comps[]) throws Exception {
		if(compList==null)
			compList = new HashMap<String, COSMOSACCompound>();

		this.comps = comps;
		this.ncomps = comps.length;

		this.VCOSMO = new double[ncomps];

		for (int i = 0; i < comps.length; i++) {
			SigmaProfileGenerator s = new SigmaProfileGenerator(SigmaProfileGenerator.FileType.GAMESS,
					this.rav, nsegments);
			
			COSMOSACCompound c2 = compList.get(comps[i].name);
			if(c2!=null){
				comps[i].charge = c2.charge;
				this.VCOSMO[i] = comps[i].Vcosmo = c2.Vcosmo;
				comps[i].area = c2.area;
				comps[i].areaMulti = c2.areaMulti;
				continue;
			}
			
			String name = comps[i].name.replace(' ','_');
			String extension = ".gout";
//			String folder = "mopac/";
			String folder = "moltest/";
			
			try {
				s.parseFile(folder + name + extension, rav);												
			} catch (FileNotFoundException e) {
				s.parseFile(folder + name.replace('-','_') + extension, rav);
			}
			
			comps[i].charge = s.getChargeDensity();
			this.VCOSMO[i] = comps[i].Vcosmo = s.getVolume();
			
			// the multi-area
			double[] sigma1 = s.getAveragedChargeDensity();

			String extensionLow = ".low" + extension;
			try {
				s.parseFile(folder + name + extensionLow, rav);
//				s.parseFile(folder + name + extension, rav*2);
			} catch (FileNotFoundException e) {
				s.parseFile(folder + name.replace('-','_') + extensionLow, rav);
//				s.parseFile(folder + name.replace('-','_') + extension, rav*2);
			}
			double[] sigma2 = s.getAveragedChargeDensity();

			// value from Klamt (COSMO-RS refinement)
			double fcorr = 0.816;
			fcorr = 0.9;

			double[] area = s.getOriginalArea();
			double[] area2 = new double[area.length];
			double[] area3 = new double[area.length];
			double[] sigmaT = new double[area.length];

			for (int m = 0; m < area.length; m++) {
				sigmaT[m] = 1000*Math.abs(sigma2[m]-fcorr*sigma1[m]);
			}
			
			// only 2 dimensions
			comps[i].areaMulti = new double[ndescriptors][];
			
			double tLimit = 0.7; // so that METHANE is entirely on the first category
			double tLimit2 = 2;
			for (int m = 0; m < area.length; m++) {
				if(sigmaT[m]>tLimit){
					if(sigmaT[m]<tLimit2)
						area2[m] = area[m];
					else
						area3[m] = area[m];
					area[m] = 0;
				}
			}
			
			
			// now a new filter for area3 and area2
			try {
				s.parseFile(folder + name + extension, rav*2);												
			} catch (FileNotFoundException e) {
				s.parseFile(folder + name.replace('-','_') + extension, rav*2);
			}
			sigma2 = s.getAveragedChargeDensity();
			for (int m = 0; m < area.length; m++) {
				double sT2 = 1000*Math.abs(sigma2[m]-fcorr*sigma1[m]);
				if(area3[m]>0 && sT2<1.5){
					area2[m]+= area3[m];
					area3[m] = 0;
				}
			}

			
			s.simpleSorting(area, sigma1);
			comps[i].areaMulti[0] = s.getSortedArea();
			s.simpleSorting(area2, sigma1);
			comps[i].areaMulti[1] = s.getSortedArea();
			s.simpleSorting(area3, sigma1);
			comps[i].areaMulti[2] = s.getSortedArea();
			
			compList.put(comps[i].name, comps[i]);
			
//			s.printProfile(System.out);
		}
		this.T = 300;

		z = new double[ncomps];
		for (int i = 0; i < ncomps; i++)
			z[i] = 1.0/ncomps;

		parametersChanged();
	}
}
