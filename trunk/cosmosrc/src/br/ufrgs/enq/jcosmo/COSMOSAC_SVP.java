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
public class COSMOSAC_SVP extends COSMOSAC {
	public String toString(){
		return "COSMO-SAC(SVP)";
	}

	public COSMOSAC_SVP(int numberOfSegments) {
		super(numberOfSegments);
		
		setSigmaHB(0.00);
//		setFpol(FPOL);
//		setCHB(0);
//		setAnorm(ANORM);
//		setCoord(COORD);
	}

	public COSMOSAC_SVP() {
		this(51);
	}

	public void setComponents(COSMOSACCompound comps[]) throws Exception {
		this.ncomps = comps.length;
		this.comps = comps;

		this.VCOSMO = new double[ncomps];
		this.area = new double[ncomps][];

		SigmaProfileGenerator s = new SigmaProfileGenerator(SigmaProfileGenerator.FileType.GAMESS_SVP,
				this.rav, nsegments);
		for (int i = 0; i < comps.length; i++) {
			String name = comps[i].name.replace(' ','_');
			String extension = ".svp.gout";
			String folder = "moltest/";
			
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

		parametersChanged();
	}
}
