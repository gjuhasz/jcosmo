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
	
		// teste COST:0.32256382240158127 NP:242
//		setAEffPrime(6.998695161135107);
////		setCoord(10.0);
//		setVnorm(91.5841666358196);
//		setAnorm(33.69882992335153);
//		setCHB(57479.50661977082);
//		setSigmaHB(0.006949893900864999);
//		setEpsilon(10.255836590044746);
		
		// teste COST:0.3725789961109445 NP:343 gamma
//		setAEffPrime(6.5);
////		setCoord(COORD);
////		setVnorm(VNORM);
//		setAnorm(44.10);
//		setCHB(19880.0);
//		setSigmaHB(0.0027);
//		setEpsilon(4.7); 
		
		// teste COST:0.3514230794222891 NP:343 lngamma%
//		setAEffPrime(6.82037037037037);
////		setCoord(COORD);
////		setVnorm(VNORM);
//		setAnorm(22.975333333333335);
//		setCHB(85580.0);
//		setSigmaHB(0.0064088888888888884);
//		setEpsilon(4.970822222222221);
		
		// teste COST:0.45866623662193606 NP:343 gamma
//		setAEffPrime(7.413137623530041);
////		setCoord(COORD);
////		setVnorm(VNORM);
//		setAnorm(31.120902318292327);
//		setCHB(65243.713860566);
//		setSigmaHB(0.0064);
//		setEpsilon(6.090955782578806);
		
// 		// optimized parameters COST:0.3455019859714253 NP:343 lngamma% am1 *
//		setAEffPrime(7.25);
//		setCoord(10.0);
//		setVnorm(66.69);
//		setAnorm(25.85);
//		setCHB(85580.0);
//		setSigmaHB(0.0064);
//		setEpsilon(4.971);

		// teste RMS=0.5536787472192798
//		setAEffPrime(6.577619440499027);//1.0496235771209508
////		setCoord(10.0);
////		setVnorm(66.69);
//		setAnorm(27.237869264957105);//3.905968470576764
//		setCHB(53921.70250431103);//20611.804905777906
//		setSigmaHB(0.0064);
//		setEpsilon(7.399713408532802);//3.715212187255123
		
		//** Article1 
//		setAEffPrime(6.5);
////		setCoord(10.0);
////		setVnorm(66.69);
//		setAnorm(27.25);
//		setCHB(53920.0);
//		setSigmaHB(0.0064);
//		setEpsilon(7.395);
		
//		//** Aarticle2 RMS=0.5536787472192798
//		setAEffPrime(7.6);
////		setCoord(10.0);
////		setVnorm(66.69);
//		setAnorm(22.33);
//		setCHB(54420.0);
//		setSigmaHB(0.0064);
////		setEpsilon(3.667);
		
		//** Aarticle3 RMS=0.6416224763232018
////		setAEffPrime(7.5);
////		setCoord(10.0);
////		setVnorm(66.69);
//		setAnorm(22.43);
//		setCHB(54778.0);
//		setSigmaHB(0.0064);
////		setEpsilon(3.667);
		
		//** Article4 RMS=0.5880605501590166
		setAEffPrime(6.2);
//		setCoord(10.0);
//		setVnorm(66.69);
		setAnorm(28.20);
		setCHB(42700.0);
		setSigmaHB(0.0064);
		setEpsilon(13.139);//+-24.359508118522637
		
	}

	public void setComponents(COSMOSACCompound comps[]) throws Exception {
		this.ncomps = comps.length;

		this.VCOSMO = new double[ncomps];
		this.sigma = new double[ncomps][];

		for (int i = 0; i < comps.length; i++) {
			SigmaProfileGenerator s = null;
			
			try {
				s = new SigmaProfileGenerator(SigmaProfileGenerator.FileType.MOPAC,
						"mopac/" + comps[i].name + ".cos");												
			} catch (Exception e) {
				s = new SigmaProfileGenerator(SigmaProfileGenerator.FileType.MOPAC,
						"mopac/" + comps[i].name.replace('-',' ') + ".cos");
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
