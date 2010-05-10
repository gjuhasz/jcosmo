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
		
		// Only non-polar elements, COST:0.08499795473312625 NP:32
//		idac/Alkane-Alkane.csv AARD:0.02157167555584826 NP:3
//		idac/Alkane-AlkylHalide.csv AARD:0.20398044102438556 NP:5
//		idac/AlkylHalide-Alkane.csv AARD:0.09639315605035846 NP:13
//		idac/Aromatic-Alkane.csv AARD:0.013181456097686173 NP:6
//		idac/CycloAlkane-AlkylHalide.csv AARD:0.06062351088595813 NP:5
		setBeta(0.6720307269285903);
		setCHB(0.0);
		setSigmaHB(0.01);
		setSigmaHB2(0.0);
		setSigmaHB3(1.0);
		setFpol(0.7508070397398836);
		setIgnoreSG(false);
		setCoord(10);
		setAnorm(79.53);
		setVnorm(66.69);
		
//		idac/Alcohol-Alkane.csv AARD:0.4982266327131054 NP:14
//		idac/Ketone-Alcohol.csv AARD:0.2816430413805225 NP:34, COST:0.3448131635730743
		setCHB(0);
		setSigmaHB(0.004);
		setCHB(1, 120000);
		
//		COST:0.27958305302063197
//		idac/Alcohol-Water.csv AARD:1.6029666923270154 NP:4
//		idac/Aldehyde-Water.csv AARD:NaN NP:0
//		idac/Alkane-Water.csv AARD:0.27997699244638746 NP:14
//		idac/Alcohol-Alkane.csv AARD:0.355200978079096 NP:4
//		idac/Alcohol-CycloAlkane.csv AARD:0.14638123370339154 NP:5
//		idac/Alkane-Alcohol.csv AARD:0.3295395129026112 NP:5
//		idac/Alkane-Alkane.csv AARD:0.019217495724087444 NP:2
//		idac/Alkane-AlkylHalide.csv AARD:0.16328428366109754 NP:5
//		idac/Alkane-Ketone.csv AARD:0.3118924529858303 NP:16
//		idac/Alkane-Phenol.csv AARD:0.23212603926156566 NP:16
//		idac/AlkylHalide-Alkane.csv AARD:0.2661814486873221 NP:10
//		idac/Aromatic-Alkane.csv AARD:0.278165824996176 NP:3
//		idac/CycloAlkane-Alcohol.csv AARD:0.10813960619617895 NP:3
//		idac/CycloAlkane-AlkylHalide.csv AARD:0.03038726189890581 NP:5
//		idac/CycloAlkane-Amine.csv AARD:1.5217135556152046 NP:1
//		idac/CycloAlkane-Phenol.csv AARD:0.1518539060387431 NP:5
//		idac/Ketone-Alkane.csv AARD:0.08374733315346959 NP:13
		setBeta(1.1423500947431098);
		setFpol(0.7726267939852054);
		setFpol(1, 0.7977790401526919);
		setCHB(0.0);
		setCHB(1, 0.3013797759750783);
		setSigmaHB(0.0025322357476481724);
		setSigmaHB2(0.0);
		setIgnoreSG(false);
		setCoord(10.0);
		setAnorm(79.53);
		setVnorm(66.69);
		
		// HB estimation, COST:0.34553669837717527 NP:24
//		idac/Alcohol-Water.csv AARD:0.52902010362079 NP:8
//		idac/Aldehyde-Water.csv AARD:0.12205295553595463 NP:2
//		idac/Alkane-Water.csv AARD:0.27261715667961706 NP:14
		setCHB(0);
		
		// nonaqueous+aqueous, COST:0.32677799151994164 NP:143
//		idac/Alcohol-Water.csv AARD:0.4988574700550016 NP:10
//		idac/Water.csv AARD:0.2687208633249313 NP:23
//		idac/Alcohol-Alkane.csv AARD:0.5953974262710329 NP:10
//		idac/Alcohol-CycloAlkane.csv AARD:0.4402448289876905 NP:10
//		idac/Alkane-Alcohol.csv AARD:0.6500254618124048 NP:9
//		idac/Alkane-Alkane.csv AARD:0.019350497916921217 NP:2
//		idac/Alkane-AlkylHalide.csv AARD:0.28924093440961784 NP:5
//		idac/Alkane-Ketone.csv AARD:0.14803814371151758 NP:16
//		idac/Alkane-Phenol.csv AARD:0.18924327367406021 NP:16
//		idac/AlkylHalide-Alkane.csv AARD:0.43796133972889306 NP:10
//		idac/Aromatic-Alkane.csv AARD:0.5537557307290863 NP:3
//		idac/CycloAlkane-Alcohol.csv AARD:0.4139295161647365 NP:6
//		idac/CycloAlkane-AlkylHalide.csv AARD:0.030019360180914084 NP:5
//		idac/CycloAlkane-Phenol.csv AARD:0.17048964584829668 NP:5
//		idac/Ketone-Alkane.csv AARD:0.2265411891359618 NP:13
		setBeta(0.7871234028286334);
		setFpol(1.5699046534819647);
		setFpol(1, 0.42939157969075425);
		setFpol(2, 0.38311458616327376);
		setSigmaHB(0.006);
		setCHB(2, 10067);
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
			String folder = "mopac/";
//			String folder = "moltest/";
			
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
			
			double tLimit = 0.3; // so that METHANE is entirely on the first category
			double tLimit2 = 1;
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
