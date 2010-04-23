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
 * COSMO-SAC (GAMAESS sigma profiles) activity model.
 *
 * @author Rafael de Pelegrini Soares
 * 
 */
public class COSMOSAC_G extends COSMOSAC {
	public String toString(){
		return "COSMO-SAC(GAMESS)";
	}

	public COSMOSAC_G(int numberOfSegments) {
		super(numberOfSegments);
		
		setBeta(1.4141655544216012);
		setCHB(21856.59194516772);
		setSigmaHB(0.0033588530289351158);
		setFpol(0.17960975053558242);
		setIgnoreSG(false);
		setCoord(19.621223483264043);
		setAnorm(77.08098812555656);
		setVnorm(66.69);
		
		setFpol(FPOL);
	}

	public COSMOSAC_G() {
		this(51);
	}

	public void setComponents(COSMOSACCompound comps[]) throws Exception {
		this.ncomps = comps.length;

		this.VCOSMO = new double[ncomps];
		this.area = new double[ncomps][];

		SigmaProfileGenerator s = new SigmaProfileGenerator(SigmaProfileGenerator.FileType.GAMESS,
				SigmaProfileGenerator.RAV, nsegments);
		for (int i = 0; i < comps.length; i++) {
			String name = comps[i].name.replace(' ','_');
			String extension = ".gout";
			String folder = "moltest/";
			
			try {
				s.parseFile(folder + name + extension);												
			} catch (Exception e) {
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
