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
 * COSMO-SAC (MOPAC sigma profiles) activity model.
 *
 * @author Rafael de Pelegrini Soares
 * 
 */
public class COSMOPACMultiAtom extends COSMOSACMulti {
	public String toString(){
		return "COSMO-SAC(MOPAC)";
	}
	
	String folder = "mopac/";
	
	public COSMOPACMultiAtom(){
		super(51, 4);
		
		// we use another averaging radius
		this.rav = COSMOSAC.RAV*1.1;
		
		// nonHB, COST:0.08398846320047755, NP:68
		folder = "mopAM1/";
		setBeta(1.9455);
		setCHB(0.0);
		setFpol(0.54815);
		setFpol(1, 0.552345);
		setFpol(2, 0.43642137);
		setFpol(2, 0.5023036);
		setIgnoreSG(false);
		setCoord(10.0);
		setAnorm(42.795026);
		setVnorm(65.027547);
		// HB part COST:0.3471202179898385 NP:81
		setSigmaHB(0.004);
		setCHB(1,2, 2935);
		setCHB(1,3, 1462);
		
		// aqueous and nonaqueous, COST:0.5459409975779058 NP:246
		folder = "mopAM1/";
		setBeta(1);
		setFpol(1.0259720937136563);
		setFpol(1, 0.8236345);
		setFpol(2, 0.9746551);
		setFpol(3, 0.8355474);
		setIgnoreSG(false);
		setCoord(10.0);
		setAnorm(42.795026);
		setVnorm(65.027547);
		setCHB(1,2, 1519);
		setCHB(1,3, 150);
		
		// aqueous and nonaqueous, vdw*1.18, COST:0.436 NP:246
		folder = "mopAM1_R1.1/";
		setBeta(1);
		setFpol(1.238143);
		setFpol(1, 0.20150);
		setFpol(2, 2.12397);
		setFpol(3, 1.10189);
		setIgnoreSG(false);
		setCoord(10.0);
		setAnorm(42.795026);
		setVnorm(65.027547);
		setCHB(1,2, 2263.);
		setCHB(1,3, 280);
	}
	
	protected void calculeDeltaW_HB(){
		for(int d=0; d<ndescriptors; ++d){
			for(int d2=0; d2<ndescriptors; ++d2){
				for(int m=0; m<nsegments; ++m){
					for(int n=0; n<nsegments; ++n){
						int ACC = n, DON = m;
						if(charge[m]>=charge[n]){
							ACC = m;
							DON = n;
						}
						// Hydrogen Bond effect:
						double hb = 0;
						if((charge[ACC] * charge[DON]) < 0)
							hb = (charge[ACC] - charge[DON]);
						hb = -hb*hb;
//						double hb = Math.max(0.0, charge[ACC] - sigmaHB)*Math.min(0.0, charge[DON] + sigmaHB);

						deltaW_HB[d][d2][m][n] = cHB[d][d2]* hb;
					}
				}
			}
		}
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
				s.parseFile(folder + name + extension);												
			} catch (FileNotFoundException e) {
				s.parseFile(folder + name.replace('-','_') + extension);
			}
			
			MolParser molParser = new MolParser();
			
			try{
				molParser.parseFile(folder + name + ".mol");
			}
			catch (FileNotFoundException e) {
				molParser.parseFile(folder + name.replace('-','_') + ".mol");
			}
			int[] bondAtom1 = molParser.getBondAtom1();
			int[] bondAtom2 = molParser.getBondAtom2();
			int[] elementType = molParser.getElementType();

			
			comps[i].charge = s.getChargeDensity();
			int[] elem = s.getElem();
			int[] atoms = s.getAtom();
			this.VCOSMO[i] = comps[i].Vcosmo = s.getVolume();
			
			// the multi-area
			double[] sigma1 = s.getAveragedChargeDensity();

			double[] area = s.getOriginalArea();

			comps[i].areaMulti = new double[ndescriptors][];
			
			// lets filter the H-[N,O,F,Cl,Br,I] atoms, the donnors
			double[] areaHDonnor = new double[area.length];
			for (int m = 0; m < area.length; m++) {
				if(elem[m]==1){
					for (int j = 0; j < bondAtom1.length; j++) {
						if(
							((bondAtom1[j]==atoms[m] && elementType[bondAtom2[j]-1]==7))
							|| ((bondAtom1[j]==atoms[m] && elementType[bondAtom2[j]-1]==8))
							|| ((bondAtom1[j]==atoms[m] && elementType[bondAtom2[j]-1]==9))
							|| ((bondAtom1[j]==atoms[m] && elementType[bondAtom2[j]-1]==17))
							|| ((bondAtom1[j]==atoms[m] && elementType[bondAtom2[j]-1]==35))
							|| ((bondAtom1[j]==atoms[m] && elementType[bondAtom2[j]-1]==53))
							){
							areaHDonnor[m] += area[m];
							area[m] = 0;
							break;
						}
						else if(
							((elementType[bondAtom1[j]-1]==7 && bondAtom2[j]==atoms[m]))
							|| ((elementType[bondAtom1[j]-1]==8 && bondAtom2[j]==atoms[m]))
							|| ((elementType[bondAtom1[j]-1]==9 && bondAtom2[j]==atoms[m]))
							|| ((elementType[bondAtom1[j]-1]==17 && bondAtom2[j]==atoms[m]))
							|| ((elementType[bondAtom1[j]-1]==35 && bondAtom2[j]==atoms[m]))
							|| ((elementType[bondAtom1[j]-1]==53 && bondAtom2[j]==atoms[m]))
							){
							areaHDonnor[m] += area[m];
							area[m] = 0;
							break;
						}
					}
				}
			}
			
			// lets filter the O-H atoms
			double[] areaOH = new double[area.length];
			for (int m = 0; m < area.length; m++) {
				if(elem[m]==8){
					for (int j = 0; j < bondAtom1.length; j++) {
						if( (bondAtom1[j]==atoms[m] && elementType[bondAtom2[j]-1]==1) ||
								(elementType[bondAtom1[j]-1]==1 && bondAtom2[j]==atoms[m])){
							areaOH[m] += area[m];
							area[m] = 0;
							break;
						}
					}
				}
			}
			
			// lets filter the N,O,F,Cl,Br,I atoms, acc
			double[] areaAcc = new double[area.length];
			for (int m = 0; m < area.length; m++) {
				if(elem[m]==7 || elem[m]==8 || elem[m]==9 || elem[m]==17 || elem[m]==35 || elem[m]==53){
					areaAcc[m] += area[m];
					area[m] = 0;
				}
			}

			s.simpleSorting(area, sigma1);
			comps[i].areaMulti[0] = s.getSortedArea();
			s.simpleSorting(areaHDonnor, sigma1);
			comps[i].areaMulti[1] = s.getSortedArea();
			s.simpleSorting(areaOH, sigma1);
			comps[i].areaMulti[2] = s.getSortedArea();
			s.simpleSorting(areaAcc, sigma1);
			comps[i].areaMulti[3] = s.getSortedArea();
			
//			s.printProfile(System.out);
		}
		this.T = 300;

		z = new double[ncomps];
		for (int i = 0; i < ncomps; i++)
			z[i] = 1.0/ncomps;

		parametersChanged();
	}
}
