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
		
		// COST:0.1974776564263834 NP:111
//		idac/Alcohol-Alkane.csv AARD:0.32655062475976615 NP:10
//		idac/Alcohol-CycloAlkane.csv AARD:0.32250821425094867 NP:15
//		idac/Alkane-Alcohol.csv AARD:0.12881474043480728 NP:9
//		idac/Alkane-Alkane.csv AARD:0.02129552686818885 NP:3
//		idac/Alkane-AlkylHalide.csv AARD:0.15866138428457155 NP:2
//		idac/Alkane-Ketone.csv AARD:0.13380342525878233 NP:12
//		idac/Alkane-Phenol.csv AARD:0.10227286077191758 NP:22
//		idac/AlkylHalide-Alkane.csv AARD:0.3711002026563018 NP:6
//		idac/Amine-Alkane.csv AARD:NaN NP:0
//		idac/Aromatic-Alkane.csv AARD:0.1565140462657024 NP:6
//		idac/CycloAlkane-Alcohol.csv AARD:0.18358899730050174 NP:6
//		idac/CycloAlkane-AlkylHalide.csv AARD:NaN NP:0
//		idac/CycloAlkane-Phenol.csv AARD:0.2351685077179408 NP:5
//		idac/Ketone-Alcohol.csv AARD:0.33495840669380617 NP:7
//		idac/Ketone-Alkane.csv AARD:0.07939862746498699 NP:8
		setBeta(1.805216580498087);
		setFpol(0.53759919833672);
		setFpol(1, 0.5240924876831323);
		setFpol(2, 0.653533925543);
		setSigmaHB(0.008);
		setSigmaHB2(0.0);
		setSigmaHB3(1.0);
		setCHB(0);
		setCHB(2, 15958);
		setIgnoreSG(false);
		setCoord(10);
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
			
			double tLimit = 1; // so that METHANE is entirely on the first category
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
			s.parseFile(folder + name + extension, rav*2);
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
