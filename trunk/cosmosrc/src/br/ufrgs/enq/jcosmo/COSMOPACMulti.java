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

import java.io.File;
import java.io.FileNotFoundException;
import java.util.HashMap;

import br.ufrgs.enq.jcosmo.SigmaProfileGenerator.FileType;



/**
 * COSMO-SAC (MOPAC sigma profiles) activity model.
 *
 * @author Rafael de Pelegrini Soares
 * 
 */
public class COSMOPACMulti extends COSMOSACMulti {
	private String folder;
	private HashMap<String, COSMOSACCompound> compList;
	private FileType type;
	private String extension;

	public String toString(){
		return "COSMO-SAC(MOPAC)";
	}
	
	public COSMOPACMulti(){
		super(21, 5);
		
		// we use another averaging radius
		this.rav = COSMOSAC.RAV;
		
		type = SigmaProfileGenerator.FileType.MOPAC;
		extension = ".cos";
		folder = "mopac/";
//		folder = "moltest/";
//		folder = "mopAM1/";
//		folder = "mopAM1c/";
//		folder = "mopRM1/";
		
//		type = SigmaProfileGenerator.FileType.GAMESS;
//		extension = ".gout";
//		folder = "gam6-31Gd/";
		
//		// article results
//		setResCorr(1);
//		setCHB(42700.7265672813);
//		setSigmaHB(0.0064);
//		setFpol(FPOL);
//		setIgnoreSG(false);
//		setAnorm(28.2);
//		setVnorm(66.69);
//		setAEff(7.5);
		
		// Estimating a base Fpol, Anorm and beta
		// nonHB.csv NP:94  COST:0.08950252780992553
		// EPS=999.0 CHARGE=0 COSWRT RSOLV=1.2 AM1 VDW(H=1.276:C=1.972:N=1.798:O=1.7632:F=1.7052:S=2.088:P=2.088:Cl=2.03:Br=2.146:I=2.2968) GNORM=0.1 RELSCF=0.1
		setBeta(1.9242298625056287);
		setCHB(0.0);
		setFpol(0.5351729083650381);
		setAnorm(52.253508045274685);
		
		setCHB(88700);
		setSigmaHB(0.009);
		setSigmaHB2(0.004);
		setBeta(2.25);
		setFpol(0.43);
		sigmaDisp = 0.0005;
		cDisp = 0;
		
		setCHB(0);
		setSigmaHB(0);
		setSigmaHB2(0);
		setCHB(1, 2, 6404);
		setCHB(1, 3, 11003);
		setAnorm(80);
		
		// All nonaqueous, organic acids removed
		// COST:0.18910604411365706, NP:309
		setBeta(1.6193867595103684);
		setFpol(0.664523609196695);
		setAnorm(157.60541396778427);
		setSigmaHB(0);
		setCHB(1, 1, 10319);
		setCHB(1, 2, 4000);
		
		// aqueous
		// COST:8259199383240307, NP:281
		setCHB(1, 4, 3051);
		setCHB(3, 1, 8681);
		setCHB(3, 2, 757);
		setCHB(3, 4, 4172);
		
		// aqueous, without some multiring aromatics
		// COST:0.6484449484937216 NP:276
		setCHB(1, 4, 3538);
		setCHB(3, 1, 8485);
		setCHB(3, 2, 372);
		setCHB(3, 4, 3973);
	}

	public void setComponents(COSMOSACCompound comps[]) throws Exception {
		this.comps = comps;
		this.ncomps = comps.length;

		this.VCOSMO = new double[ncomps];
		
		if(compList==null)
			compList = new HashMap<String, COSMOSACCompound>();

		SigmaProfileGenerator s = new SigmaProfileGenerator(type,
				this.rav, nsegments);
		for (int i = 0; i < comps.length; i++) {
			
			String name = comps[i].name.replace(' ','_');
			
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
			int[] atoms = s.getAtom();

			double[] area0 = s.getOriginalArea();
			double[] area1 = new double[area0.length];
			double[] area2 = new double[area0.length];
			double[] area3 = new double[area0.length];
			double[] area4 = new double[area0.length];

			comps[i].areaMulti = new double[ndescriptors][];
			
			MolParser molParser = new MolParser();
			File molFile = new File(folder + name + ".mol");
			if(!molFile.exists())
				molFile = new File(folder + name + ".mol+1");
			if(!molFile.exists())
				molFile = new File(folder + name + ".mol-1");

			molParser.parseFile(molFile.getPath());
			
			// H2O in its own group
			if(comps[i].name.equals("WATER")){
				for (int m = 0; m < area0.length; m++) {
					if(sigma1[m]>0 && molParser.matchType(atoms[m], 8, 0)){
						area4[m] += area0[m];
						area0[m] = 0;
					}
					else
						if(sigma1[m]<0 && molParser.matchType(atoms[m], 1, 0)){
						area3[m] += area0[m];
						area0[m] = 0;
					}
				}
			}

			// lets filter the H-[N,O,...] atoms (HB-donor)
			for (int m = 0; m < area0.length; m++) {
				int boundedType[] = {7, 8, 9, 17, 35, 53};
				if(sigma1[m]<0 && molParser.matchType(atoms[m], 1, boundedType)){
					area1[m] += area0[m];
					area0[m] = 0;
				}
			}
			
			// lets filter the [N,O,...]-H atoms (HB-acceptor bounded to H)
			for (int m = 0; m < area0.length; m++) {
				int atomType[] = {7, 8, 9, 17, 35, 53};
				if(sigma1[m]>0 && molParser.matchType(atoms[m], atomType, 1)){
					area1[m] += area0[m];
					area0[m] = 0;
				}
			}
			
			// lets filter the [N, O, ...] atoms
			for (int m = 0; m < area0.length; m++) {
				int atomType[] = {7, 8, 9, 17, 35, 53};
//				All on group2
//				if(sigma1[m]>0 && molParser.matchType(atoms[m], atomType, 0)){
//					area2[m] += area0[m];
//					area0[m] = 0;
//				}

				// Oxygen goes to area2 only if double bounded
				if(sigma1[m]>0 && molParser.matchBondType(atoms[m], atomType, 2)){
					area2[m] += area0[m];
					area0[m] = 0;
				}
				if(sigma1[m]>0 && molParser.matchBondType(atoms[m], atomType, 3)){
					area2[m] += area0[m];
					area0[m] = 0;
				}
				// single bounded Oxygen go to area1
				// FIXME: detect automatically ACETATE and OXIDE oxygens
				if(!comps[i].name.endsWith("ATE") && !comps[i].name.endsWith("OXIDE")){
					if(sigma1[m]>0 && molParser.matchType(atoms[m], 8, 0)){
						area1[m] += area0[m];
						area0[m] = 0;
					}
				}
				if(sigma1[m]>0 && molParser.matchType(atoms[m], atomType, 0)){
					area2[m] += area0[m];
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
			s.simpleSorting(area4, sigma1);
			comps[i].areaMulti[4] = s.getSortedArea();
			
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
