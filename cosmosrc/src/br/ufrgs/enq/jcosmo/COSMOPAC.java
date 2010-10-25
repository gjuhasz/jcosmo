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
	
	String folder = "mopac/";

	public COSMOPAC() {
		super(51);
		
		// we use another averaging radius
		this.rav = COSMOSAC.RAV;
		
		folder = "mopAM1/";
		folder = "mopRM1/";
		folder = "mopPOA1/";
		
//		// article results
		setBeta(1);
		setCHB(42700.7265672813);
		setSigmaHB(0.0064);
		setFpol(FPOL);
		setIgnoreSG(false);
		setAnorm(28.2);
		setVnorm(66.69);
		
		setAnorm(ANORM);
		setBeta(1.29150929);
		setFpol(0.43043818);
//		sigmaDisp = 9.999989892401997E-4;
//		cDisp = 22.3570;
//		
//		// Only non-HB and low polarizability (RAV*1.1),  COST:0.08096689998854072
////		idac/Alkane-Alkane.csv AARD:0.07719095226428561 NP:7
////		idac/Alkane-AlkylHalide.csv AARD:0.21358093090737124 NP:5
////		idac/Alkane-Ketone.csv AARD:0.11663264840429993 NP:21
////		idac/AlkylHalide-Alkane.csv AARD:0.09508975499887927 NP:35
////		idac/Aromatic-Alkane.csv AARD:0.041880708942964404 NP:46
////		idac/CycloAlkane-AlkylHalide.csv AARD:0.06457926168719777 NP:5
//		setBeta(1.6513070950256865);
////		setCHB(0.0);
////		setSigmaHB(0.006);
//		setSigmaHB2(0.0);
//		setSigmaHB3(1.0);
//		setFpol(0.6900883503832824);
//		setIgnoreSG(false);
//		setCoord(10.0);
//		setAnorm(51.404961223433276);
//		setVnorm(66.69);
//		
//		
//		// now adjusting the HB terms with all nonaqueous
////		setCHB(9.392658327026367E8);
////		setSigmaHB(0.020400000000000012);
//		
//		// nonHB, COST:0.10516900155174044 NP:165
//		setBeta(1.9871344869851588);
//		setCHB(0.0);
//		setSigmaHB(0.01);
//		setSigmaHB2(0.0);
//		setSigmaHB3(1.0);
//		setFpol(0.5360550192460967);
//		setIgnoreSG(false);
//		setCoord(10.0);
//		setAnorm(46.78140757190759);
//		setVnorm(66.69);
//		
//		// nonHB, RM1, COST:0.112, NP:68
//		setBeta(2.144806645259847);
//		setCHB(0.0);
//		setSigmaHB(2.4440796188024025);
//		setSigmaHB2(0.0);
//		setSigmaHB3(1.0);
//		setFpol(0.5049604559424442);
//		setIgnoreSG(false);
//		setCoord(10.0);
//		setAnorm(81.05566549258984);
//		setVnorm(66.69);
//		
//		// eletrostatic
//		setSigmaHB(0.008);
//		setCHB(0);
//		setFpol(151);
		
		
//		// COST:0.661595084764779
//		this.rav = COSMOSAC.RAV*0.6;
//		setBeta(1.0833791891812146);
//		setCHB(30723.010871759267);
//		setSigmaHB(0.002433385217034843);
//		setSigmaHB2(0.0);
//		setSigmaHB3(1.0);
//		setFpol(0.08819078172495727);
//		setIgnoreSG(false);
//		setCoord(10.0);
//		setAnorm(32.5848146849939);
//		setVnorm(66.69);
//
//		// COST:0.6187205854133636
//		this.rav = COSMOSAC.RAV*0.01;
//		setBeta(1.0282683136600688);
//		setCHB(28462.21891389604);
//		setSigmaHB(0.002433386633495701);
//		setSigmaHB2(0.0);
//		setSigmaHB3(1.0);
//		setFpol(0.07184423274107313);
//		setIgnoreSG(false);
//		setCoord(10.0);
//		setAnorm(37.76131078724174);
//		setVnorm(66.69);
//		
//		// :H=1.3:C=2.0:N=1.83:O=1.79:F=1.72:S=2.16:P=2.12:Cl=2.07:Br=2.18:I=2.32
//		// COST:0.559
//		folder = "mopAM1c/";
//		this.rav = COSMOSAC.RAV*0.01;
//		setBeta(1.2612819680772671);
//		setCHB(27744.33994865579);
//		setSigmaHB(0.002585296448832151);
//		setSigmaHB2(0.0);
//		setSigmaHB3(1.0);
//		setFpol(0.07847475358309021);
//		setIgnoreSG(false);
//		setCoord(10.0);
//		setAnorm(34.69401869906048);
//		setVnorm(66.69);
//		
//		// EPS=999.0 CHARGE=0 COSWRT RSOLV=1.2 AM1 VDW(:H=1.44:C=2.21:N=2.015:O=1.976:F=1.911:S=2.34:P=2.34:Cl=2.275:Br=2.405:I=2.574) GNORM=0.1 RELSCF=0.1
//		// COST:0.6023343630822742
//		this.rav = COSMOSAC.RAV*0.01;
//		setBeta(0.7736216294236089);
//		setCHB(34835.9827671406);
//		setSigmaHB(0.001195408597551723);
//		setSigmaHB2(0.0);
//		setSigmaHB3(1.0);
//		setFpol(0.04342506402019329);
//		setIgnoreSG(false);
//		setCoord(10.0);
//		setAnorm(55.93073585292235);
//		setVnorm(66.69);
//		
//		// EPS=999.0 CHARGE=0 COSWRT RSOLV=2.4 AM1 VDW(:H=1.44:C=2.21:N=2.015:O=1.976:F=1.911:S=2.34:P=2.34:Cl=2.275:Br=2.405:I=2.574) GNORM=0.1 RELSCF=0.1
//		// COST:0.5679029405844803
//		setBeta(1.182348955348154);
//		setCHB(49271.95312893527);
//		setSigmaHB(0.0020234303032980686);
//		setSigmaHB2(0.0);
//		setSigmaHB3(1.0);
//		setFpol(0.09037784897659495);
//		setIgnoreSG(false);
//		setCoord(10.0);
//		setAnorm(27.417788010046856);
//		setVnorm(66.69);
//		
//		// EPS=5 CHARGE=0 COSWRT RSOLV=1.2 AM1 VDW(:H=1.44:C=2.21:N=2.015:O=1.976:F=1.911:S=2.34:P=2.34:Cl=2.275:Br=2.405:I=2.574) GNORM=0.1 RELSCF=0.1
//		// COST:0.5454579234422399
//		setBeta(1.1685916291195768);
//		setCHB(97657.79992765447);
//		setSigmaHB(0.0013826219894970337);
//		setSigmaHB2(0.0);
//		setSigmaHB3(1.0);
//		setFpol(0.15711050803756574);
//		setIgnoreSG(false);
//		setCoord(10.0);
//		setAnorm(30.281177263334282);
//		setVnorm(66.69);
//		
////		// EPS=5 CHARGE=0 COSWRT RSOLV=1.2 AM1 VDW(:H=1.3:C=2.0:N=1.83:O=1.79:F=1.72:S=2.16:P=2.12:Cl=2.07:Br=2.18:I=2.32) GNORM=0.1 RELSCF=0.1
////		// COST:0.546049465502303
////		setBeta(1.0565472966585658);
////		setCHB(59652.74717744731);
////		setSigmaHB(0.0014648860401811998);
////		setSigmaHB2(0.0);
////		setSigmaHB3(1.0);
////		setFpol(0.11609313112102097);
////		setIgnoreSG(false);
////		setCoord(10.0);
////		setAnorm(23.659474255246465);
////		setVnorm(66.69);
////		
//////		EPS=5 CHARGE=0 COSWRT RSOLV=2.2 AM1 VDW(:H=1.3:C=2.0:N=1.83:O=1.79:F=1.72:S=2.16:P=2.12:Cl=2.07:Br=2.18:I=2.32) GNORM=0.1 RELSCF=0.1
//////		COST:0.6391708893913468
////		setBeta(0.8331366330666732);
////		setCHB(66843.29958814124);
////		setSigmaHB(0.0011012137936171128);
////		setSigmaHB2(0.0);
////		setSigmaHB3(1.0);
////		setFpol(0.16764606120299697);
////		setIgnoreSG(false);
////		setCoord(10.0);
////		setAnorm(43.769203298804996);
////		setVnorm(66.69);
//		
//		// EPS=5 CHARGE=0 COSWRT RSOLV=0.6 AM1 VDW(:H=1.3:C=2.0:N=1.83:O=1.79:F=1.72:S=2.16:P=2.12:Cl=2.07:Br=2.18:I=2.32) GNORM=0.1 RELSCF=0.1
//		// COST:0.535390567158958
//		setBeta(0.9495328456872831);
//		setCHB(59622.270498422586);
//		setSigmaHB(8.713395874295679E-4);
//		setSigmaHB2(0.0);
//		setSigmaHB3(1.0);
//		setFpol(0.09208773507305704);
//		setIgnoreSG(false);
//		setCoord(10.0);
//		setAnorm(20.69697773316389);
//		setVnorm(66.69);
		
		// nonHB COST:0.1228614057740516 NP:196
		// RSOLV=1.2 RM1 VDW(H=1.276:C=1.972:N=1.898:O=1.8632:F=1.7052:S=2.088:P=2.088:Cl=2.43:Br=2.146:I=2.2968)
		folder = "profiles/RM1/";
		rav = 1.15*RAV;
		setBeta(1.4174886027020652);
		setFpol(0.7258152987045381);
		setCHB(0);
		setIgnoreSG(false);
		setCoord(10.0);
		setAnorm(53.96283430888596);
		setRPower(0.8305611275321202);
		setVnorm(86.01343230125369);
		
		// nonHB COST:0.1295918375389806 NP:196
		// RSOLV=1.2 AM1 EXTERNAL=POA1.rm1 VDW(H=1.392:C=1.972:N=1.798:O=1.7632:F=1.7052:S=2.088:P=2.088:Cl=2.03:Br=2.146:I=2.2968)
		folder = "profiles/POA1/";
		rav = 1.15*RAV;
		setBeta(1.7112762508999033);
		setCHB(0.0);
		setFpol(0.6884343905835382);
		setIgnoreSG(false);
		setCoord(10.0);
		setAnorm(56.502077074611016);
		setRPower(0.8152912289882516);
		setVnorm(86.01343230125369);
		
		// Non energetic study
		// idac/Athermal.csv AARD:0.03634304855190082 NP:146
		folder = "profiles/RM1/";
		setAnorm(225.53381478405623);
		setRPower(0.6229842090038551);
		// idac/Athermal.csv AARD:0.036278854768369694 NP:146
		folder = "profiles/RM1_1.18/";
		setAnorm(124.18305909611365);
		setRPower(0.6372590453541258);
		// idac/nonHB.csv AARD:0.1271605266268788 NP:269
		setBeta(1.399460298613211);
		setFpol(0.73377065315896);
		
		folder = "profiles/RM1/";
		
		// idac/nonHB.csv AARD:0.10439292160859014 NP:269
		rav = 1.30*RAV;
		setBeta(1.3205203534195487);
		setFpol(0.843371077585408);
		setAnorm(64.61423251936094);
		setRPower(0.7700935807819056);
	}

	public void setComponents(COSMOSACCompound comps[]) throws Exception {
		this.comps = comps;
		this.ncomps = comps.length;

		this.VCOSMO = new double[ncomps];
		this.area = new double[ncomps][];

		SigmaProfileGenerator s = new SigmaProfileGenerator(SigmaProfileGenerator.FileType.MOPAC,
				this.rav, nsegments);
		for (int i = 0; i < comps.length; i++) {
			
			String name = comps[i].name.replace(' ','_');
			String extension = ".cos";
			
			try {
				s.parseFile(folder + name + extension);												
			} catch (Exception e) {
				s.parseFile(folder + name.replace('-','_') + extension);
			}
			
			comps[i].charge = s.getChargeDensity();
			this.VCOSMO[i] = comps[i].Vcosmo = s.getVolume();
			this.area[i] = comps[i].area = s.getSortedArea();
			
//			s.printProfile(System.out);
		}
		super.setComponents(comps);
	}
}
