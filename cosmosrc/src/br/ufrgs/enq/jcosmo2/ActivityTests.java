
package br.ufrgs.enq.jcosmo2;


public class ActivityTests {

	public static void main(String[] args) throws Exception {
		
		Compound comps[] = new Compound[2];
		comps[0] = new Compound();
		comps[1] = new Compound();

		comps[0].name = "BENZENE";
		comps[1].name = "ETHYLENE CARBONATE";
		
//		comps[0].name = "CARBON TETRACHLORIDE";
//		comps[1].name = "N-HEXANE";
		
//		comps[0].name = "N,N-DIMETHYLFORMAMIDE";
//		comps[1].name = "N-HEXANE";

//		comps[0].name = "N-OCTANE";
//		comps[1].name = "ETHYL ACETATE";
		
//		comps[0].name = "CYCLOHEXYLAMINE";
//		comps[1].name = "N-OCTANE";
		
		comps[0].name = "ACETONE";
		comps[1].name = "N-HEPTANE";

//		comps[0].name = "N-PENTANE";
//		comps[1].name = "ETHYL ACETATE";
//
//		comps[0].name = "N-PENTANE";
//		comps[1].name = "METHYL ETHYL KETONE";
		
//		comps[0].name = "TOLUENE";
//		comps[1].name = "ACETONITRILE";
		
		COSMOSAC2 cosmosac = new COSMOSAC2();
		cosmosac.setComponents(comps);

		cosmosac.setTemperature(298.15);

		int n = 21;
		double x[] = new double[2];
		double lnGamma[] = new double[2];

		System.out.println("X1\tGamma1\tGamma2\tlnGamma1\tlnGamma2");

		for (int i = 0; i < n; i++) {
			x[0] = (double)i/(n-1);
			x[1] = 1-x[0];
			
			cosmosac.setComposition(x);
			cosmosac.activityCoefficient(lnGamma);
			
			System.out.println(x[0] +"\t"+ Math.exp(lnGamma[0]) +"\t"+ Math.exp(lnGamma[1]) +"\t"+ lnGamma[0] +"\t"+ lnGamma[1]);
		}
	}
}
