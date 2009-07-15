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
		
		// optimization based on VLE experiments
//		setAEffPrime(2.032);
//		setCoord(0.4954);
//		setVnorm(9.11500);
//		setAnorm(32.71439);
//		setCHB(31747.393);
//		setSigmaHB(0.006491);
		
		// teste COST:0.33605919189470757 NP:264
//		setAEffPrime(7.104872197041376);
//		setCoord(10.0);
//		setVnorm(50.356158728594586);
//		setAnorm(40.34455758285999);
//		setCHB(38909.18398032317);
//		setSigmaHB(0.007610431584682716);
		
		// teste  COST:0.3349320427103625 NP:264
//		setAEffPrime(6.987770682010926);
//		setCoord(10.831896133137302);
//		setVnorm(56.70929194847275);
//		setAnorm(41.20603421625156);
//		setCHB(39338.90939677377);
//		setSigmaHB(0.007689283637684755);
		
		// teste  COST:0.36202518040747816 NP:310 
//		setAEffPrime(7.5);
//		setCoord(10.0);
//		setVnorm(39.312074722512);
//		setAnorm(44.47306992385012);
//		setCHB(33385.50979424656);
//		setSigmaHB(0.0072033250851981315);
		
		// teste COST:0.5048779355889912 NP:285
//		setAEffPrime(7.5);
//		setCoord(10.0);
//		setVnorm(42.07135995186812);
//		setAnorm(44.5202111100251);
//		setCHB(30559.688907926902);
//		setSigmaHB(0.007010968843097385);
		
//		 teste 0.51
//		setAEffPrime(7.5);
//		setCoord(10.0);
//		setVnorm(29.059646484375016);
//		setAnorm(36.52438549804689);
//		setCHB(52973.5185546875);
//		setSigmaHB(0.008460908203124996);
		
		// teste COST:0.464473140176675
//		setAEffPrime(7.5);
//		setCoord(10.0);
		setVnorm(38.15387540947493);
		setAnorm(40.685560114024895);
		setCHB(29662.893349642523);
		setSigmaHB(0.0066128361885796965);
		
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
