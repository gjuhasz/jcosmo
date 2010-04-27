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
public class COSMOPACMulti extends COSMOSACMulti {
	public String toString(){
		return "COSMO-SAC(MOPAC)";
	}
	
	public COSMOPACMulti(){
		super(51, 3);
		
		// we use another averaging radius
		this.rav = COSMOSAC.RAV*1.1;
		
//		// article results
//		setResCorr(1);
//		setCHB(42700.7265672813);
//		setSigmaHB(0.0064);
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
		setCHB(0.0);
		setSigmaHB(0.006);
		setSigmaHB2(0.0);
		setSigmaHB3(1.0);
		setFpol(0.6900883503832824);
		setIgnoreSG(false);
		setCoord(10.0);
		setAnorm(51.404961223433276);
		setVnorm(66.69);
		
		
		// now adjusting the HB with 3 descriptors, all nonaqueous, COST:0.4591648062037473 NP:325
		setBeta(1.6513070950256865);
		setBeta(1, 0.9698277112144933);
		setBeta(2, 2.5773245598710877);
		setCHB(132965.77713782096);
		setCHB(1, 109115.48920332415);
		setCHB(2, 78041.74607782175);
		setSigmaHB(0.0055);
		setSigmaHB2(0.0);
		setSigmaHB3(1.0);
		setFpol(0.6900883503832824);
		setFpol(1, 0.6269643704378696);
		setFpol(2, 0.8938514253981169);
		setIgnoreSG(false);
		setCoord(10.0);
		setAnorm(51.404961223433276);
		setVnorm(66.69);
		
		setSigmaHB(0.0055);
		setCHB(5e8);
		
		
		setBeta(1.6513070950256865);
		setBeta(1, 2.5966623022422692);
		setBeta(2, 3.186655778007093);
		setCHB(0.0);
		setCHB(1, 0.0);
		setCHB(2, 0.0);
		setSigmaHB(0.0055);
		setSigmaHB2(0.0);
		setSigmaHB3(1.0);
		setFpol(0.6900883503832824);
		setFpol(1, 0.4430225813238091);
		setFpol(2, 1.0318244680427964);
		setIgnoreSG(false);
		setCoord(10.0);
		setAnorm(51.404961223433276);
		setVnorm(66.69);
		
		
//		idac/Alcohol-Alkane.csv AARD:0.4982266327131054 NP:14
//		idac/Ketone-Alcohol.csv AARD:0.2816430413805225 NP:34, COST:0.3448131635730743
		setBeta(1.6513070950256865);
		setBeta(1, 4.783547453257033);
		setBeta(2, 0.7103895052380395);
		setCHB(4.190793926262297E7);
		setCHB(1, 4.190793926262297E7);
		setCHB(2, 4.190793926262297E7);
		setSigmaHB(0.0055);
		setSigmaHB2(0.0);
		setSigmaHB3(1.0);
		setFpol(0.6900883503832824);
		setFpol(1, 0.8244733229051346);
		setFpol(2, 3.342358909084158);
		setIgnoreSG(false);
		setCoord(10.0);
		setAnorm(51.404961223433276);
		setVnorm(66.69);
	}

	public void setComponents(COSMOSACCompound comps[]) throws Exception {
		this.comps = comps;
		this.ncomps = comps.length;

		this.VCOSMO = new double[ncomps];

		SigmaProfileGenerator s = new SigmaProfileGenerator(SigmaProfileGenerator.FileType.MOPAC,
				this.rav, nsegments);
		for (int i = 0; i < comps.length; i++) {
			
			String name = comps[i].name.replace(' ','_');
			String extension = ".cos";
			
			try {
				s.parseFile("mopac/" + name + extension);												
			} catch (Exception e) {
				s.parseFile("mopac/" + name.replace('-','_') + extension);
			}
			
			comps[i].charge = s.getChargeDensity();
			this.VCOSMO[i] = comps[i].Vcosmo = s.getVolume();
			
			// the multi-area
			double[] sigma1 = s.getAveragedChargeDensity();

			s.averageCharges(rav*2);
			double[] sigma2 = s.getAveragedChargeDensity();

			// value from Klamt (COSMO-RS refinement)
			double fcorr = 0.816;

			double[] area = s.getOriginalArea();
			double[] area2 = new double[area.length];
			double[] area3 = new double[area.length];
			double[] sigmaT = new double[area.length];

			s.simpleSorting(area, sigma1);

			for (int m = 0; m < area.length; m++) {
				//					sigmaT[m] = 1000*(sigma2[m]-fcorr*sigma1[m]);
				sigmaT[m] = 1000*Math.abs(sigma2[m]-fcorr*sigma1[m]);
			}
			
			// only 2 dimensions
			comps[i].areaMulti = new double[ndescriptors][];
			
			double tLimit = 2;
			double tLimit2 = 3;
			for (int m = 0; m < area.length; m++) {
				if(sigmaT[m]>tLimit){
					if(sigmaT[m]>tLimit2)
						area3[m] = area[m];
					else
						area2[m] = area[m];
					area[m] = 0;
				}
			}
			
			s.simpleSorting(area, sigma1);
			comps[i].areaMulti[0] = s.getSortedArea();
			s.simpleSorting(area2, sigma1);
			comps[i].areaMulti[1] = s.getSortedArea();
			s.simpleSorting(area3, sigma1);
			comps[i].areaMulti[2] = s.getSortedArea();
			
//			s.printProfile(System.out);
		}
		this.T = 300;

		z = new double[ncomps];
		for (int i = 0; i < ncomps; i++)
			z[i] = 1.0/ncomps;

		parametersChanged();
	}
	
	protected void calculeDeltaW_HB(){
		for (int d = 0; d < ndescriptors; d++) {
			for(int m=0; m<nsegments; ++m){
				for(int n=0; n<nsegments; ++n){
					int ACC = n, DON = m;
					if(charge[m]>=charge[n]){
						ACC = m;
						DON = n;
					}
					// Hydrogen Bond effect:
					double hb = Math.max(0.0, charge[ACC] - sigmaHB)*Math.min(0.0, charge[DON] + sigmaHB);

					// Klamt, Fluid Phase Equilib. 2000
					// double cHBT_c = 1.5;
					double cHBT = 1; // Math.max(0, 1 + cHBT_c * (298.15/T - 1));

					hb = -hb*hb;

					deltaW_HB[d][m][n] = cHB[d]*cHBT* hb;
				}
			}
		}
	}
}
