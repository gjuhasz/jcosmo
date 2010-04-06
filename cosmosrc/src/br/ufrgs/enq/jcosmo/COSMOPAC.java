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

import br.ufrgs.enq.jcosmo.SigmaProfileGenerator.FileType;


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
		// article results
		setResCorr(1);
		setCHB(42700.7265672813);
		setSigmaHB(0.0064);
		setFpol(FPOL);
		setIgnoreSG(false);
		setAnorm(28.2);
		setVnorm(66.69);		
		
		// IDAC with all systems but amines, COST:
//		setResCorr(1);
//		setCHB(CHB);
//		setSigmaHB(sigmaHB);
//		setFpol(FPOL);
//		setIgnoreSG(false);
//		setAnorm(ANORM);
//		setVnorm(VNORM);
	}

	public void setComponents(COSMOSACCompound comps[]) throws Exception {
		this.comps = comps;		
		this.ncomps = comps.length;

		this.VCOSMO = new double[ncomps];
		this.sigma = new double[ncomps][];

		SigmaProfileGenerator s = new SigmaProfileGenerator(SigmaProfileGenerator.FileType.MOPAC);
		for (int i = 0; i < comps.length; i++) {
			
			String name = comps[i].name.replace(' ','_');
			String extension = ".cos";
			
			try {
				s.parseFile("mopac/" + name + extension);												
			} catch (Exception e) {
				s.parseFile("mopac/" + name.replace('-','_') + extension);
			}
			
			this.charge = comps[i].charge = s.getChargeDensity();
			this.VCOSMO[i] = comps[i].Vcosmo = s.getVolume();
			this.sigma[i] = comps[i].sigma = s.getSigmaProfile();
		}
		this.compseg = charge.length;
		this.T = 300;

		z = new double[ncomps];
		for (int i = 0; i < ncomps; i++)
			z[i] = 1.0/ncomps;

		parametersChanged();
	}
}
