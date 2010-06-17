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


/**
 * COSMO-SAC (GAMAESS sigma profiles) activity model.
 *
 * @author Rafael de Pelegrini Soares
 * 
 */
public class COSMOSAC_G extends COSMOSAC {
	private String folder = "mopac/";

	public String toString(){
		return "COSMO-SAC(GAMESS)";
	}

	public COSMOSAC_G(int numberOfSegments) {
		super(numberOfSegments);
		
//		folder = "moltest/";
//		folder = "gam6-31+G2d,p/";
		folder = "gam6-31Gd/";
//		folder = "gamSTO3/";
		
		this.rav = RAV;

		
		setFpol(FPOL);
		setCHB(0);
		setAnorm(ANORM);
		setCoord(COORD);

//		Only non-HB elements, COST:0.0952035855018203 NP:50
//		idac/Alkane-Alkane.csv AARD:0.027880685952315475 NP:3
//		idac/Alkane-AlkylHalide.csv AARD:0.17044823778183177 NP:5
//		idac/Alkane-Ketone.csv AARD:0.06669337926724342 NP:18
//		idac/AlkylHalide-Alkane.csv AARD:0.17058206591526598 NP:13
//		idac/Aromatic-Alkane.csv AARD:0.018372007393053225 NP:6
//		idac/CycloAlkane-AlkylHalide.csv AARD:0.059203260051550174 NP:5
		setBeta(4.056604339857187);
		setCHB(0.0);
		setSigmaHB(0.01);
		setSigmaHB2(0.0);
		setSigmaHB3(1.0);
		setFpol(0.29880411235834803);
		setIgnoreSG(false);
		setCoord(7.2);
		setAnorm(211.90488545095045);
		setVnorm(66.69);
		
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
		
		setBeta(2.716723022971033);
		setCHB(0.0);
		setSigmaHB(0.01);
		setSigmaHB2(0.0);
		setSigmaHB3(1.0);
		setFpol(1.7074753121035136);
		setIgnoreSG(false);
		setCoord(10.0);
		setAnorm(ANORM);
		setVnorm(66.69);
		
		
		setFpol(0.52);
		setCHB(8000);
		setBeta(0.97);
		setSigmaHB(0.001);
		setSigmaHB2(0.0088);
	}

	public COSMOSAC_G() {
		this(11);
	}

	public void setComponents(COSMOSACCompound comps[]) throws Exception {
		this.ncomps = comps.length;
		this.comps = comps;

		this.VCOSMO = new double[ncomps];
		this.area = new double[ncomps][];

		SigmaProfileGenerator s = new SigmaProfileGenerator(SigmaProfileGenerator.FileType.GAMESS,
				this.rav, nsegments);
		for (int i = 0; i < comps.length; i++) {
			String name = comps[i].name.replace(' ','_');
			String extension = ".gout";
			
			try {
				s.parseFile(folder + name + extension);												
			} catch (FileNotFoundException e) {
				s.parseFile(folder + name.replace('-','_') + extension);
			}
			
			comps[i].charge = s.getChargeDensity();
			this.VCOSMO[i] = comps[i].Vcosmo = s.getVolume();
			this.area[i] = comps[i].area = s.getSortedArea();
		}
		this.T = 300;

		z = new double[ncomps];
		for (int i = 0; i < ncomps; i++)
			z[i] = 1.0/ncomps;

		super.setComponents(comps);
		parametersChanged();
	}
}
