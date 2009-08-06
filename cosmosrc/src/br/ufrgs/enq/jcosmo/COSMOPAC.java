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
		
		// teste COST:0.37800887481053136 NP:242
//		setAEffPrime(6.406312590256911);
////		setCoord(10.0);
//		setVnorm(49.08102528394544);
//		setAnorm(48.8073386838568);
//		setCHB(24991.230313953296);
//		setSigmaHB(0.0022658291813892955);
//		setEpsilon(3.2756158636530275);
		
		// teste COST:0.32256382240158127 NP:242
//		setAEffPrime(6.998695161135107);
////		setCoord(10.0);
//		setVnorm(91.5841666358196);
//		setAnorm(33.69882992335153);
//		setCHB(57479.50661977082);
//		setSigmaHB(0.006949893900864999);
//		setEpsilon(10.255836590044746);
		
		// teste COST:0.33753272552129526 NP:242
//		setAEffPrime(7.092291285998105);
////		setCoord(10.0);
//		setVnorm(21.501990883549087);
//		setAnorm(32.252260681169446);
//		setCHB(42020.95995384622);
//		setSigmaHB(0.00886888028547375);
//		setEpsilon(3.8698365749213774);
//		setSigmaGaussian(true);
		
		// teste COST:0.2887708953007761 NP:242
//		setAEffPrime(6.472620447898134);
//		setCoord(10.0);
//		setVnorm(53.159972679322344);
//		setAnorm(48.679351183766364);
//		setCHB(22145.117295919925);
//		setSigmaHB(0.0021335704762264537);
//		setEpsilon(2.99944948223479);
				
		// teste COST:0.2563595580149685 NP:242
//		setAEffPrime(5.833182929456752);
////		setCoord(COORD);
////		setVnorm(VNORM);
////		setAnorm(ANORM);
//		setCHB(19049.08512470434);
//		setSigmaHB(0.0014294669608294712);
//		setEpsilon(2.812848033405693);
		
		// teste COST:0.25779779158503163 NP:242
//		setAEffPrime(6.01151508164739);
////		setCoord(COORD);
////		setVnorm(VNORM);
////		setAnorm(ANORM);
//		setCHB(17014.10714339448);
//		setSigmaHB(0.0014599629489520833);
////		setEpsilon(EPSILON);
		
		// teste COST:0.2568374014591271 NP:242
////		setAEffPrime(AEFFPRIME);
////		setCoord(COORD);
////		setVnorm(VNORM);
////		setAnorm(ANORM);
//		setCHB(16967.584194174266);
//		setSigmaHB(0.0014936105325412134);
////		setEpsilon(EPSILON);
		
		// teste COST:0.4448006278656969 NP:343
//		setAEffPrime(8.109259259259261);
//		setCoord(13.555555555555554);
////		setVnorm(VNORM);
//		setAnorm(37.114000000000004);
//		setCHB(131222.66666666666);
////		setSigmaHB(SIGMAHB);
//		setEpsilon(6.491948148148148);
		
		// teste COST:0.3725789961109445 NP:343 gamma
//		setAEffPrime(6.5);
////		setCoord(COORD);
////		setVnorm(VNORM);
//		setAnorm(44.10);
//		setCHB(19880.0);
//		setSigmaHB(0.0027);
//		setEpsilon(4.7); 
		
		// teste COST:0.47677915807855176 NP:343 lngamma
//		setAEffPrime(6.5);
////		setCoord(COORD);
////		setVnorm(VNORM);
//		setAnorm(34.50);
//		setCHB(15590.5);
//		setSigmaHB(0.0027);
//		setEpsilon(10.25);
		
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
		
		// optimized parameters COST:0.3455019859714253 NP:343 lngamma%
		setAEffPrime(7.25);
		setCoord(10.0);
		setVnorm(66.69);
		setAnorm(25.85);
		setCHB(85580.0);
		setSigmaHB(0.0064);
		setEpsilon(4.971);
		
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
