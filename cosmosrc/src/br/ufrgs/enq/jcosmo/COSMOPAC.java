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
		// TODO: needs a reparametrization based on experiments (VLE, LLE, etc)!!
		
		// optimization based on IDAC (opt)
//		setAEffPrime(7.18165706857084);
//		setCoord(13.3026);
//		setVnorm(82.629);
//		setAnorm(55.78477);
//		setCHB(69460.875);
//		setSigmaHB(0.008250);

		// optimization based on IDAC (Direct)
//		setAEffPrime(5.08);
//		setCoord(17.43);
//		setVnorm(66.69);
//		setAnorm(53.87);
//		setCHB(35594.636);
//		setSigmaHB(0.01464);
		
////		// optimization based on VLE experiments
//		setAEffPrime(2.032);
//		setCoord(0.4954);
//		setVnorm(9.11500);
//		setAnorm(32.71439);
//		setCHB(31747.393);
//		setSigmaHB(0.006491);
	}

	public void setComponents(COSMOSACCompound comps[]) throws Exception {
		this.ncomps = comps.length;

		this.VCOSMO = new double[ncomps];
		this.sigma = new double[ncomps][];

		for (int i = 0; i < comps.length; i++) {
			SigmaProfileGenerator s = new SigmaProfileGenerator(SigmaProfileGenerator.FileType.MOPAC,
					"mopac/" + comps[i].name + ".cos");
			
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
