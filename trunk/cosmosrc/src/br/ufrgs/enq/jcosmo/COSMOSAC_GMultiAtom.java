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
public class COSMOSAC_GMultiAtom extends COSMOSACMulti {
	public String toString(){
		return "COSMO-SAC(MOPAC)";
	}
	
	String folder = "mopac/";
	
	public COSMOSAC_GMultiAtom(){
		super(51, 4);
		

		// nonHB, COST:0.14029708507644462 NP:68
		folder = "gamSTO3/";
		setBeta(3.3849341216002857);
		setCHB(0);
		setCHB(1, 2, 3153.28125);
		setCHB(1, 3, 4000.28125);
		setFpol(1.2386513272071946);
		setFpol(1, 1.2024559836598478);
		setFpol(2, 1.5822106093430);
		setFpol(3, 1.46507954035040);
		setIgnoreSG(false);
		setCoord(10.0);
		setAnorm(46.78140757190759);
		setVnorm(66.69);
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
		
		SigmaProfileGenerator s = new SigmaProfileGenerator(SigmaProfileGenerator.FileType.GAMESS,
				this.rav, nsegments);
		for (int i = 0; i < comps.length; i++) {
			
			String name = comps[i].name.replace(' ','_');
			String extension = ".gout";
			
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
