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
 * COSMO-SAC (MOPAC sigma profiles) activity model.
 *
 * @author Rafael de Pelegrini Soares
 * 
 */
public class COSMOPACMulti extends COSMOSACMulti {
	private String folder;
	private HashMap<String, COSMOSACCompound> compList;

	public String toString(){
		return "COSMO-SAC(MOPAC)";
	}
	
	public COSMOPACMulti(){
		super(11, 4);
		
		// we use another averaging radius
		this.rav = COSMOSAC.RAV*1.1;
		
//		folder = "mopac/";
//		folder = "moltest/";
		folder = "mopAM1/";
//		folder = "mopAM1c/";
//		folder = "mopRM1/";
		
//		// article results
//		setResCorr(1);
//		setCHB(42700.7265672813);
//		setSigmaHB(0.0064);
//		setFpol(FPOL);
//		setIgnoreSG(false);
//		setAnorm(28.2);
//		setVnorm(66.69);
//		setAEff(7.5);
		
//		// nonHB.csv (RAV*1.1), NP:68  COST:0.07479337893069214
//		setBeta(1.87);
//		setCHB(0.0);
//		setFpol(0.4345827604481827);
//		setFpol(1, 0.6182737464687393);
//		setFpol(2, 0.4948283183309885);
//		setFpol(3, 0.6838034897132979);
//		setIgnoreSG(false);
//		setCoord(10.0);
//		setAnorm(51.404961223433276);
//		setVnorm(66.69);
//		
//		// EPS=999.0 CHARGE=0 COSWRT RSOLV=1.2 AM1 VDW(:H=1.3:C=2.0:N=1.83:O=1.72:F=1.72:S=2.16:P=2.12:Cl=2.07:Br=2.18:I=2.32) GNORM=0.1 RELSCF=0.1
//		// idac/nonaqueous.csv AARD:0.4176077942821443 NP:166
//		// idac/aqueous.csv AARD:0.5244229674120215 NP:80
//		setSigmaHB(0.008388239);
//		setSigmaHB2(0);
//		setCHB(1, 2, 797);
//		setCHB(2, 1, -1.74);
//		setCHB(1, 3, 19300.);
//		setCHB(3, 1, 19401.);
//		
//		// non-aqueous only, COST:0.3137
//		setBeta(1.7502858561732015);
//		setFpol(0.3457623637926366);
//		setFpol(1, 0.6935853680360622);
//		setFpol(2, 0.4941785425065802);
//		setFpol(3, 0.7683764418119046);
//		setAnorm(60.819943053930515);
//		setSigmaHB(0.0088);
//		setSigmaHB2(0.0088);
//		setCHB(1, 2, 100);
//		setCHB(2, 1, 100);
//		setCHB(1, 3, 107.);
//		setCHB(3, 1, 100.);
		
//		// Groups: others; H; O,N,etc; O-H, COST:0.4532669707828741
//		// idac/nonaqueous.csv AARD:0.440472175632788 NP:166
//		// idac/aqueous.csv AARD:0.48003946329413677 NP:80
//		setCHB(0);
//		setBeta(2.3069996666555173);
////		setFpol(0.47323656400537706);
////		setFpol(1, 0.5592631058376497);
////		setFpol(2, 0.4208543674397245);
////		setFpol(3, 1.2156411479251445);
//		setAnorm(112.85141665207878);
//		
//		setFpol(0.47323656400537706);
//		setFpol(1,1, 0.5592631058376497);
//		setFpol(2,2, 0.4208543674397245);
//		setFpol(3,3, 1.2156411479251445);
//		for (int i = 0; i < 4; i++) {
//			for (int j = 0; j < 4; j++) {
//				setFpol(i, j, Math.sqrt(getFpol(i, i)*getFpol(j, j)));
//			}
//		}
		
		
		// Groups: others; H; O,N,etc; O-H, COST:0.3648489587564063 NP:246
		// idac/nonaqueous.csv AARD:0.3361917496296972 NP:166
		// idac/aqueous.csv AARD:0.42431266769432785 NP:80
		setCHB(0);
		setBeta(1.824074284569755);
		setAnorm(149.125);
		
		// 1.8286285937091067
		// 0.21303723739603353 0.47895356566983804 0.2941563710679771 1.3628976010855547
		// 0.38646287577471794 0.7987742630400895 0.6463922165674696 0.8134684979374178
		// 0.3888527719462525 0.5066232026537906 0.2058210361459763 0.85719378375602
		// 0.9547030809105133 0.9840153879560437 1.0559812655063003 1.5948680633435508
		// 149.52779145885367  COST:0.3648489587564063 NP:246

		setFpol(0,0, 0.21303723);
		setFpol(0,1, 0.47895356566983804);
		setFpol(0,2, 0.2941563710679771);
		setFpol(0,3, 1.3628976010855547);
		
		setFpol(1,0, 0.38646287577471794);
		setFpol(1,1, 0.7987742630400895);
		setFpol(1,2, 0.6463922165674696);
		setFpol(1,3, 0.8134684979374178);
		
		setFpol(2,0, 0.3888527719462525);
		setFpol(2,1, 0.5066232026537906);
		setFpol(2,2, 0.2058210361459763);
		setFpol(2,3, 0.85719378375602);

		setFpol(3,0, 0.9547030809105133);
		setFpol(3,1, 0.9840153879560437);
		setFpol(3,2, 1.0559812655063003);
		setFpol(3,3, 1.5948680633435508);
	}

	public void setComponents(COSMOSACCompound comps[]) throws Exception {
		this.comps = comps;
		this.ncomps = comps.length;

		this.VCOSMO = new double[ncomps];
		
		if(compList==null)
			compList = new HashMap<String, COSMOSACCompound>();

		SigmaProfileGenerator s = new SigmaProfileGenerator(SigmaProfileGenerator.FileType.MOPAC,
				this.rav, nsegments);
		for (int i = 0; i < comps.length; i++) {
			
			String name = comps[i].name.replace(' ','_');
			String extension = ".cos";
			
			COSMOSACCompound c2 = compList.get(comps[i].name);
			if(c2!=null){
				comps[i].charge = c2.charge;
				this.VCOSMO[i] = comps[i].Vcosmo = c2.Vcosmo;
				comps[i].area = c2.area;
				comps[i].areaMulti = c2.areaMulti;
				continue;
			}
			
			try {
				s.parseFile(folder + name + extension);												
			} catch (FileNotFoundException e) {
				s.parseFile(folder + name.replace('-','_') + extension);
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
