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
		
		// COST:0.13071578444923296 NP:92
//		idac/Alkane-Alcohol.csv AARD:0.1457625504704864 NP:9
//		idac/Alkane-Alkane.csv AARD:0.021306165795884433 NP:3
//		idac/Alkane-AlkylHalide.csv AARD:0.12951704985191428 NP:5
//		idac/Alkane-Ketone.csv AARD:0.10124178607952884 NP:18
//		idac/Alkane-Phenol.csv AARD:0.10800138084182385 NP:22
//		idac/AlkylHalide-Alkane.csv AARD:0.17793718759396865 NP:13
//		idac/Aromatic-Alkane.csv AARD:0.17738558436957882 NP:6
//		idac/CycloAlkane-Alcohol.csv AARD:0.18484509202466623 NP:6
//		idac/CycloAlkane-AlkylHalide.csv AARD:0.035588755861995794 NP:5
//		idac/CycloAlkane-Phenol.csv AARD:0.22791833281982857 NP:5
		setBeta(1.8677387992263745);
		setFpol(0.58049990);
		setFpol(1, 0.292276224);
		setFpol(2, 0.421617194);
		setSigmaHB(0.01);
		setSigmaHB2(0.0);
		setSigmaHB3(1.0);
		setCHB(0);
		setCHB(2, 113220);
		setIgnoreSG(false);
		setCoord(10);
		setAnorm(107.16496305);
		setVnorm(66.69);
		
		// COST:0.1678591113051289 NP:131
//		idac/Alcohol-Alkane.csv AARD:0.3494957698383428 NP:10
//		idac/Alcohol-CycloAlkane.csv AARD:0.35208720003494387 NP:15
//		idac/Alkane-Alcohol.csv AARD:0.1439600000593901 NP:9
//		idac/Alkane-Alkane.csv AARD:0.021291524439811763 NP:3
//		idac/Alkane-AlkylHalide.csv AARD:0.12670897616856733 NP:5
//		idac/Alkane-Ketone.csv AARD:0.10665831131060521 NP:18
//		idac/Alkane-Phenol.csv AARD:0.10531834779726368 NP:22
//		idac/AlkylHalide-Alkane.csv AARD:0.18794664139966455 NP:13
//		idac/Aromatic-Alkane.csv AARD:0.20539215138547348 NP:6
//		idac/CycloAlkane-Alcohol.csv AARD:0.14719374417565048 NP:6
//		idac/CycloAlkane-AlkylHalide.csv AARD:0.030938558548510752 NP:5
//		idac/CycloAlkane-Phenol.csv AARD:0.23228294906313335 NP:5
//		idac/Ketone-Alkane.csv AARD:0.07920339041020255 NP:14
		setBeta(1.6827599663147546);
		setFpol(0.6407712814407316);
		setFpol(1, 0.278454604684);
		setFpol(2, 0.46127354229949);
		setSigmaHB(0.01);
		setSigmaHB2(0.0);
		setSigmaHB3(1.0);
		setCHB(0);
		setCHB(2, 113220);
		setIgnoreSG(false);
		setCoord(10);
		setAnorm(66196.627);
		setVnorm(66.69);
		
		// COST:0.21093663696082796 NP:145
//		idac/Alcohol-Alkane.csv AARD:0.3413052989309958 NP:10
//		idac/Alcohol-CycloAlkane.csv AARD:0.3036033029873327 NP:15
//		idac/Alkane-Alcohol.csv AARD:0.20581962679001597 NP:9
//		idac/Alkane-Alkane.csv AARD:0.0317541110073166 NP:3
//		idac/Alkane-AlkylHalide.csv AARD:0.14190790532239403 NP:5
//		idac/Alkane-Ketone.csv AARD:0.11258850685166771 NP:18
//		idac/Alkane-Phenol.csv AARD:0.11121564605853264 NP:22
//		idac/AlkylHalide-Alkane.csv AARD:0.1437458698552522 NP:13
//		idac/Aromatic-Alkane.csv AARD:0.18032107248279508 NP:6
//		idac/CycloAlkane-Alcohol.csv AARD:0.11871154660128531 NP:6
//		idac/CycloAlkane-AlkylHalide.csv AARD:0.0286507435446411 NP:5
//		idac/CycloAlkane-Phenol.csv AARD:0.21159016223902838 NP:5
//		idac/Ketone-Alcohol.csv AARD:0.6522044711099412 NP:14
//		idac/Ketone-Alkane.csv AARD:0.10727904180436301 NP:14
		setBeta(1.675128096644348);
		setFpol(0.658228278044782);
		setFpol(1, 0.10625869937826941);
		setFpol(2, 0.4001722977102460);
		setSigmaHB(0.01);
		setSigmaHB2(0.0);
		setSigmaHB3(1.0);
		setCHB(0);
		setCHB(2, 32935.717);
		setIgnoreSG(false);
		setCoord(10);
		setAnorm(66196.627);
		setVnorm(66.69);
	}
	
	protected void calculeDeltaW_HB(){
		for(int d=0; d<ndescriptors; ++d){
			for(int d2=0; d2<ndescriptors; ++d2){
				for(int m=0; m<nsegments; ++m){
					for(int n=0; n<nsegments; ++n){
						int ACC = n, DON = m;
						if(charge[m]>=charge[n]){
							ACC = m;
							DON = n;
						}
						// Hydrogen Bond effect:
						double hb = Math.min(0.0, charge[ACC]*charge[DON]);

						// Klamt, Fluid Phase Equilib. 2000
						// double cHBT_c = 1.5;
						double cHBT = 1; // Math.max(0, 1 + cHBT_c * (298.15/T - 1));
						
						hb = -Math.pow(Math.abs(hb), 2)*1000;

						deltaW_HB[d][d2][m][n] = Math.sqrt(cHB[d]*cHB[d2])*cHBT* hb;
					}
				}
			}
		}
	}


	public void setComponents(COSMOSACCompound comps[]) throws Exception {
		if(compList==null)
			compList = new HashMap<String, COSMOSACCompound>();

		this.comps = comps;
		this.ncomps = comps.length;

		this.VCOSMO = new double[ncomps];

		SigmaProfileGenerator s = new SigmaProfileGenerator(SigmaProfileGenerator.FileType.GAMESS,
				this.rav, nsegments);
		for (int i = 0; i < comps.length; i++) {
			
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
				s.parseFile(folder + name + extension);												
			} catch (FileNotFoundException e) {
				s.parseFile(folder + name.replace('-','_') + extension);
			}
			
			comps[i].charge = s.getChargeDensity();
			this.VCOSMO[i] = comps[i].Vcosmo = s.getVolume();
			
			// the multi-area
			double[] sigma1 = s.getAveragedChargeDensity();

			extension = ".low" + extension;
			try {
				s.parseFile(folder + name + extension);												
			} catch (FileNotFoundException e) {
				s.parseFile(folder + name.replace('-','_') + extension);
			}
			double[] sigma2 = s.getAveragedChargeDensity();

			// value from Klamt (COSMO-RS refinement)
			double fcorr = 0.816;
			fcorr = 0.9;

			double[] area = s.getOriginalArea();
			double[] area2 = new double[area.length];
			double[] area3 = new double[area.length];
			double[] sigmaT = new double[area.length];

			s.simpleSorting(area, sigma1);

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
