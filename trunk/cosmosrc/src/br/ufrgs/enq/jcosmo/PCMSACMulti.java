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
 * COSMO-SAC (MOPAC sigma profiles) activity model.
 *
 * @author Rafael de Pelegrini Soares
 * 
 */
public class PCMSACMulti extends COSMOSACMulti {
	// Store already loaded compounds
	static HashMap<String, COSMOSACCompound> compList;  

	public String toString(){
		return "COSMO-SAC(MOPAC)";
	}
	
	public PCMSACMulti(){
		super(51, 3);
		
		this.rav = RAV;
		
		// Results for single-domain without HB
		setBeta(1.5313645197328083);
		setCHB(0.0);
		setFpol(0.7323758674505356);
		setIgnoreSG(false);
		setCoord(10.0);
		setAnorm(361.45773761532666);
		setVnorm(66.69);
		
		// results with only 
		setBeta(1.5313645197328083);
//		setBeta(1, 4.5258695349368985);
//		setBeta(2, 0.11097188871453148);
		setFpol(0.7323758674505356);
//		setFpol(1, 0.3829759018403384);
//		setFpol(2, 564.152422356673);
		setIgnoreSG(false);
		setCoord(10.0);
		setAnorm(361.45773761532666);
		setVnorm(66.69);
		
//		setBeta(0, 2.151684613190274);
//		setBeta(1, 2.846675352269272);
//		setBeta(2, 1.093915115292504);
//		setFpol(0, 0.6603110781125799);
//		setFpol(1, 1.0575437973491533);
//		setFpol(2, 1.001611743451837);
		
		setBeta(1.6522265674626846);
		setBeta(1, 1.6522265674626846);
		setBeta(2, 1.6522265674626846);
		setCHB(36.97335297330511);
//		setCHB(1, 11275.446145976268);
//		setCHB(2, 258171.92085283238);
		setSigmaHB(0.0085);
		setSigmaHB2(0.0);
		setSigmaHB3(1.0);
		setFpol(0.6757162680146256);
//		setFpol(1, 0.02456359490938178);
//		setFpol(2, 0.6734395049672226);
		setIgnoreSG(false);
		setCoord(10.0);
		setAnorm(79.53);
		setVnorm(66.69);
	}

	public void setComponents(COSMOSACCompound comps[]) throws Exception {
		if(compList==null)
			compList = new HashMap<String, COSMOSACCompound>();

		this.comps = comps;
		this.ncomps = comps.length;

		this.VCOSMO = new double[ncomps];

		SigmaProfileGenerator s = new SigmaProfileGenerator(SigmaProfileGenerator.FileType.GAMESS_PCM,
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
			String extension = ".pcm.gout";
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

			s.averageCharges(rav*2);
			double[] sigma2 = s.getAveragedChargeDensity();

			// value from Klamt (COSMO-RS refinement)
			double fcorr = 0.816;

			double[] area = s.getOriginalArea();
			double[] area2 = new double[area.length];
			double[] area3 = new double[area.length];
			double[] sigmaT = new double[area.length];

			s.simpleSorting(area, sigma1);

			for (int m = 0; m < area.length; m++) {
				//					sigmaT[m] = 1000*(sigma2[m]-fcorr*sigma1[m]);
				sigmaT[m] = 1000*Math.abs(sigma2[m]-fcorr*sigma1[m]);
			}
			
			// only 2 dimensions
			comps[i].areaMulti = new double[ndescriptors][];
			
			double tLimit = 1;
			double tLimit2 = 2;
			for (int m = 0; m < area.length; m++) {
				if(sigmaT[m]>tLimit){
					if(sigmaT[m]>tLimit2)
						area3[m] = area[m];
					else
						area2[m] = area[m];
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
