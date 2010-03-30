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
	
	public COSMOSAC_G() {
		// TODO: needs a reparametrization based on experiments (VLE, LLE, etc)!!
	}

	public void setComponents(COSMOSACCompound comps[]) throws Exception {
		this.ncomps = comps.length;

		this.VCOSMO = new double[ncomps];
		this.sigma = new double[ncomps][];

		for (int i = 0; i < comps.length; i++) {
			SigmaProfileGenerator s = null;
			String name = comps[i].name.replace(' ','_');
			
			try {
				s = new SigmaProfileGenerator(SigmaProfileGenerator.FileType.GAMESS,
						"mopac/" + name + ".gout");												
			} catch (Exception e) {
				s = new SigmaProfileGenerator(SigmaProfileGenerator.FileType.GAMESS,
						"mopac/" + name.replace('-','_') + ".gout");
			}
			
			this.charge = s.getChargeDensity();
			this.VCOSMO[i] = s.getVolume();
			this.sigma[i] = s.getSigmaProfile();
		}
		this.compseg = charge.length;
		this.T = 300;

		z = new double[ncomps];
		for (int i = 0; i < ncomps; i++)
			z[i] = 1.0/ncomps;

		parametersChanged();
	}
}
