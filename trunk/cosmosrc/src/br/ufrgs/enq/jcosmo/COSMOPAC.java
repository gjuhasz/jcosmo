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



/**
 * COSMO-SAC (MOPAC sigma profiles) activity model.
 *
 * @author Rafael de Pelegrini Soares
 * 
 */
public class COSMOPAC extends COSMOSAC {
	public String toString(){
		return "COSMO-SAC(MOPAC)";
	}
	
	public COSMOPAC() {
		
		// we use another averaging radius
		this.rav = COSMOSAC.RAV*1.1;
		
//		// article results
//		setResCorr(1);
		setCHB(42700.7265672813);
		setSigmaHB(0.0064);
//		setFpol(FPOL);
//		setIgnoreSG(false);
//		setAnorm(28.2);
//		setVnorm(66.69);
//		setAEff(7.5);
		
		// Only non-HB and low polarizability (RAV*1.1),  COST:0.08096689998854072
//		idac/Alkane-Alkane.csv AARD:0.07719095226428561 NP:7
//		idac/Alkane-AlkylHalide.csv AARD:0.21358093090737124 NP:5
//		idac/Alkane-Ketone.csv AARD:0.11663264840429993 NP:21
//		idac/AlkylHalide-Alkane.csv AARD:0.09508975499887927 NP:35
//		idac/Aromatic-Alkane.csv AARD:0.041880708942964404 NP:46
//		idac/CycloAlkane-AlkylHalide.csv AARD:0.06457926168719777 NP:5
		setBeta(1.6513070950256865);
//		setCHB(0.0);
//		setSigmaHB(0.006);
		setSigmaHB2(0.0);
		setSigmaHB3(1.0);
		setFpol(0.6900883503832824);
		setIgnoreSG(false);
		setCoord(10.0);
		setAnorm(51.404961223433276);
		setVnorm(66.69);
		
		
		// now adjusting the HB terms with all nonaqueous
//		setCHB(9.392658327026367E8);
//		setSigmaHB(0.020400000000000012);
		
		// nonHB, COST:0.10516900155174044 NP:165
		setBeta(1.9871344869851588);
		setCHB(0.0);
		setSigmaHB(0.01);
		setSigmaHB2(0.0);
		setSigmaHB3(1.0);
		setFpol(0.5360550192460967);
		setIgnoreSG(false);
		setCoord(10.0);
		setAnorm(46.78140757190759);
		setVnorm(66.69);

	}

	public void setComponents(COSMOSACCompound comps[]) throws Exception {
		this.comps = comps;
		this.ncomps = comps.length;

		this.VCOSMO = new double[ncomps];
		this.area = new double[ncomps][];

		SigmaProfileGenerator s = new SigmaProfileGenerator(SigmaProfileGenerator.FileType.MOPAC,
				this.rav, nsegments);
		for (int i = 0; i < comps.length; i++) {
			
			String name = comps[i].name.replace(' ','_');
			String extension = ".cos";
			String folder = "moltest/";
			
			try {
				s.parseFile(folder + name + extension);												
			} catch (Exception e) {
				s.parseFile(folder + name.replace('-','_') + extension);
			}
			
			comps[i].charge = s.getChargeDensity();
			this.VCOSMO[i] = comps[i].Vcosmo = s.getVolume();
			this.area[i] = comps[i].area = s.getSortedArea();
			
//			s.printProfile(System.out);
		}
		this.T = 300;

		z = new double[ncomps];
		for (int i = 0; i < ncomps; i++)
			z[i] = 1.0/ncomps;

		parametersChanged();
	}
}
