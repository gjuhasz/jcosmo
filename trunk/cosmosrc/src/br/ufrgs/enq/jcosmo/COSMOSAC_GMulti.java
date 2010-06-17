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
import java.util.HashMap;



/**
 * COSMO-SAC (GAMESS sigma profiles) activity model.
 *
 * @author Rafael de Pelegrini Soares
 * 
 */
public class COSMOSAC_GMulti extends COSMOSACMulti {
	private HashMap<String, COSMOSACCompound> compList;
	private String folder = "mopac/";
	
	public String toString(){
		return "COSMO-SAC(GAMESS)";
	}
	
	public COSMOSAC_GMulti(){
		super(51, 4);
		
//		this.rav = 0.75;
		
//		folder = "moltest/";
//		folder = "gam6-31+G2d,p/";
		folder = "gam6-31Gd/";
//		folder = "gamSTO3/";
		
//		// NONAQUEOUS, COST:0.2736488469786023 NP:167
////		idac/Alcohol-Alkane.csv AARD:0.3239770780689189 NP:10
////		idac/Alcohol-CycloAlkane.csv AARD:0.31659960715621416 NP:15
////		idac/Alkane-Alcohol.csv AARD:0.12249019700050272 NP:9
////		idac/Alkane-Alkane.csv AARD:0.021363199212819123 NP:3
////		idac/Alkane-AlkylHalide.csv AARD:0.28308602663247073 NP:5
////		idac/Alkane-Amine.csv AARD:0.31007765605690885 NP:4
////		idac/Alkane-CarboxilicAcid.csv AARD:NaN NP:0
////		idac/Alkane-Ketone.csv AARD:0.2026584124991165 NP:18
////		idac/Alkane-Phenol.csv AARD:0.13640569891573637 NP:22
////		idac/Alkene-Amine.csv AARD:0.22498699148647958 NP:3
////		idac/AlkylHalide-Alkane.csv AARD:0.3535428023306032 NP:13
////		idac/Amine-Alkane.csv AARD:0.6165306940215931 NP:5
////		idac/Aromatic-Alkane.csv AARD:0.1232432229575623 NP:6
////		idac/CycloAlkane-Alcohol.csv AARD:0.1432363420945877 NP:6
////		idac/CycloAlkane-AlkylHalide.csv AARD:0.029828597244639533 NP:5
////		idac/CycloAlkane-Amine.csv AARD:0.18564250391234083 NP:1
////		idac/CycloAlkane-CarboxilicAcid.csv AARD:0.7034224039183667 NP:2
////		idac/CycloAlkane-Phenol.csv AARD:0.14393537997518346 NP:5
////		idac/CarboxilicAcid-Alkane.csv AARD:2.558912369738963 NP:2
////		idac/CarboxilicAcid-CycloAlkane.csv AARD:2.128247009122318 NP:2
////		idac/CycloAlkene-Amine.csv AARD:0.010337077996353564 NP:3
////		idac/Ketone-Alcohol.csv AARD:0.2906249588540736 NP:14
////		idac/Ketone-Alkane.csv AARD:0.09654562496685164 NP:14
//		setBeta(1.704710195);
//		setFpol(0.49456570084204576);
//		setFpol(1, 0.11917051421257);
//		setFpol(2, 0.57751444534);
//		setSigmaHB(0.0048793);
//		setSigmaHB2(0.0);
//		setSigmaHB3(1.0);
//		setCHB(0);
////		setCHB(2, 12318.276);
//		setIgnoreSG(false);
//		setCoord(10);
//		setAnorm(107.16496305);
//		setVnorm(66.69);
//		
//		// COST:0.44667417604299436 NP:206
//		setBeta(0.8234057258502678);
//		setBeta(1, 0.8234057258502678);
//		setBeta(2, 0.8234057258502678);
//		setCHB(0.0);
////		setCHB(1, 0.0);
////		setCHB(2, 11658.463520632726);
//		setSigmaHB(0.0080);
//		setSigmaHB2(0.0);
//		setFpol(1.590397183396032);
//		setFpol(1, 0.26171380795455323);
//		setFpol(2, 0.43386856060593726);
//		setIgnoreSG(false);
//		setCoord(10.0);
//		setAnorm(107.16496305);
//		setVnorm(66.69);
//		
//		setBeta(1.4446959148112937);
//		setBeta(1, 1.4446959148112937);
//		setBeta(2, 1.4446959148112937);
//		setCHB(40000.0);
//		setCHB(0, 0, 40000.0);
//		setCHB(0, 1, 40000.0);
//		setCHB(0, 2, 40000.0);
//		setCHB(1, 1, 40000.0);
//		setCHB(1, 2, 40000.0);
//		setCHB(2, 2, 40000.0);
//		setSigmaHB(0.1);
//		setSigmaHB2(0.0);
//		setFpol(0.851319604859719);
//		setFpol(1, 0.05729418931280367);
//		setFpol(2, 0.1932782622569627);
//		setIgnoreSG(false);
//		setCoord(10.0);
//		setAnorm(53.21629458731964);
//		setVnorm(66.69);
//		
//		
//		setBeta(1.717360793290768);
//		setSigmaHB(0.00);
//		setCHB(0);
//		setFpol(0.4815843058165584);
//		setFpol(1, 0.7145865251377206);
//		setFpol(2, 0.5481833851742162);
//		setFpol(3, 1.092353511);
////		setCHB(0, 1, 3622);
////		setCHB(1, 0, 4600);
////		setCHB(1, 2, 1000);
////		setCHB(2, 1, 1000);
////		setCHB(1, 3, 2700);
////		setCHB(3, 1, 2700);
	}
	
	public void setComponents(COSMOSACCompound comps[]) throws Exception {
		if(compList==null)
			compList = new HashMap<String, COSMOSACCompound>();

		this.comps = comps;
		this.ncomps = comps.length;

		this.VCOSMO = new double[ncomps];

		for (int i = 0; i < comps.length; i++) {
			SigmaProfileGenerator s = new SigmaProfileGenerator(SigmaProfileGenerator.FileType.GAMESS,
					this.rav, nsegments);
			
			COSMOSACCompound c2 = compList.get(comps[i].name);
			if(c2!=null){
				comps[i].charge = c2.charge;
				this.VCOSMO[i] = comps[i].Vcosmo = c2.Vcosmo;
				comps[i].area = c2.area;
				comps[i].areaMulti = c2.areaMulti;
				continue;
			}
			
			String name = comps[i].name.replace(' ','_');
			String extension = ".gout";
			
			try {
				s.parseFile(folder + name + extension, rav);												
			} catch (FileNotFoundException e) {
				name = name.replace('-','_');
				s.parseFile(folder + name + extension, rav);
			}
			
			comps[i].charge = s.getChargeDensity();
			this.VCOSMO[i] = comps[i].Vcosmo = s.getVolume();
			
			// the multi-area
			double[] sigma1 = s.getAveragedChargeDensity();
			int[] elem = s.getElem();
			int[] atoms = s.getAtom();

			double[] area0 = s.getOriginalArea();
			double[] area1 = new double[area0.length];
			double[] area2 = new double[area0.length];
			double[] area3 = new double[area0.length];

			comps[i].areaMulti = new double[ndescriptors][];
			
			MolParser molParser = new MolParser();
			molParser.parseFile(folder + name + ".mol");
			int[] bondAtom1 = molParser.getBondAtom1();
			int[] bondAtom2 = molParser.getBondAtom2();
			int[] elementType = molParser.getElementType();

			// lets filter the O-H atoms
			for (int m = 0; m < area0.length; m++) {
				if(elem[m]==8){
					for (int j = 0; j < bondAtom1.length; j++) {
						if( (bondAtom1[j]==atoms[m] && elementType[bondAtom2[j]-1]==1) ||
								(elementType[bondAtom1[j]-1]==1 && bondAtom2[j]==atoms[m])){
							area3[m] += area0[m];
							area0[m] = 0;
							break;
						}
					}
				}
			}
			
			for (int m = 0; m < area0.length; m++) {
				if(elem[m]==1){
					area1[m] = area0[m];
					area0[m] = 0;
				}
				else if( (elem[m]==7 || elem[m]==8 || elem[m]==17)){ //  && sigma1[m] > 0){
					area2[m] = area0[m];
					area0[m] = 0;
				}
			}
			
			s.simpleSorting(area0, sigma1);
			comps[i].areaMulti[0] = s.getSortedArea();
			s.simpleSorting(area1, sigma1);
			comps[i].areaMulti[1] = s.getSortedArea();
			s.simpleSorting(area2, sigma1);
			comps[i].areaMulti[2] = s.getSortedArea();
			s.simpleSorting(area3, sigma1);
			comps[i].areaMulti[3] = s.getSortedArea();
			
			compList.put(comps[i].name, comps[i]);
			
//			s.printProfile(System.out);
		}
		this.T = 300;

		z = new double[ncomps];
		for (int i = 0; i < ncomps; i++)
			z[i] = 1.0/ncomps;

		parametersChanged();
	}
}
