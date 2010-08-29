package br.ufrgs.enq.jcosmo.test;

import br.ufrgs.enq.jcosmo.COSMOSAC;
import br.ufrgs.enq.jcosmo.COSMOSACCompound;
import br.ufrgs.enq.jcosmo.SigmaProfileGenerator;

public class DebyeHuckelTest {
	public static void main(String[] args) throws Exception {
		
		double ms = 36.46;
		
		SigmaProfileGenerator s = new SigmaProfileGenerator(SigmaProfileGenerator.FileType.MOPAC);
		String folder = "mopac/";
		String extension = ".cos";
		
		COSMOSACCompound water = new COSMOSACCompound();
		water.name = "WATER";
		
		s.parseFile(folder + "WATER" + extension);
		water.charge = s.getChargeDensity();
		water.Vcosmo = s.getVolume();
		water.area = s.getSortedArea();
		
		COSMOSACCompound hcl = new COSMOSACCompound();
		hcl.name = "[H3O+1][CL-1]";

		s.parseFile(folder + "H3O+1" + extension);
		hcl.charge = s.getChargeDensity();
		hcl.Vcosmo = s.getVolume();
		hcl.area = s.getSortedArea();
		
		s.parseFile(folder + "CL-1" + extension);
		hcl.Vcosmo+= s.getVolume();
		for (int i = 0; i < hcl.area.length; i++) {
			hcl.area[i] += s.getSortedArea()[i];
		}
		
		COSMOSACCompound comps[] = {hcl, water};
		COSMOSAC cosmosac = new COSMOSAC();
		cosmosac.setBeta(1.1449668413958016);
		cosmosac.setCHB(47138.69690169665);
		cosmosac.setCHB(0);
		cosmosac.setSigmaHB(0.007351579606927401);
		cosmosac.setSigmaHB2(0.007351579606927401);
		cosmosac.setSigmaHB3(1.0);
		cosmosac.setFpol(1.5991543095613814);
		cosmosac.setIgnoreSG(false);
		cosmosac.setCoord(10.0);
		cosmosac.setAnorm(79.53);
		cosmosac.setVnorm(66.69);
		
		cosmosac.setComponents(comps);
		
		double lnGamma[] = new double[2];
		double z[] = {0, 1};
		cosmosac.setComposition(z);
		
		cosmosac.activityCoefficient(lnGamma);
		cosmosac.activityCoefficient(lnGamma);
		double lnGammaInf = lnGamma[0];
		
		int N = 100;
		System.out.println("Mi\tlnGamma+-");
		for (int i = 0; i <= N; i++) {
			double Mi = 9.0*i/N;
			
			z[0] = Mi*ms/1000;
			z[1] = 1-z[0];
			
			cosmosac.setComposition(z);
			cosmosac.activityCoefficient(lnGamma);
			
			System.out.println(Mi + "\t" + (lnGamma[0] - lnGammaInf));
		}
	}
}
