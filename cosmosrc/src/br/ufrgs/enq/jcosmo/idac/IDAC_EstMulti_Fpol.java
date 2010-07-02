package br.ufrgs.enq.jcosmo.idac;

import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math.optimization.ConvergenceChecker;
import org.apache.commons.math.optimization.CostException;
import org.apache.commons.math.optimization.CostFunction;
import org.apache.commons.math.optimization.DirectSearchOptimizer;
import org.apache.commons.math.optimization.NelderMead;
import org.apache.commons.math.optimization.PointCostPair;

import br.ufrgs.enq.direct.ObjectiveFunction;
import br.ufrgs.enq.jcosmo.COSMOPACMulti;
import br.ufrgs.enq.jcosmo.COSMOSACMulti;

public class IDAC_EstMulti_Fpol implements CostFunction, ObjectiveFunction {

	List<IDACExperimentsMulti> experiments;
	static int ndesc;

	public IDAC_EstMulti_Fpol() throws Exception {
		String modelClass = COSMOPACMulti.class.getName();
//		String modelClass = COSMOPACMultiAtom.class.getName();
//		String modelClass = PCMSACMulti.class.getName();
//		String modelClass = COSMOSAC_GMulti.class.getName();
//		String modelClass = COSMOSAC_GMultiAtom.class.getName();
		
		ndesc = 5;
		
		experiments = new ArrayList<IDACExperimentsMulti>();

//		experiments.add(new IDACExperimentsMulti("idac/Alcohol-Water.csv", modelClass));
//		experiments.add(new IDACExperimentsMulti("idac/Aldehyde-Water.csv", modelClass));
//		experiments.add(new IDACExperimentsMulti("idac/Alkane-Water.csv", modelClass));
////		experiments.add(new IDACExperimentsMulti("idac/Alkene-Water.csv", modelClass));
////		experiments.add(new IDACExperimentsMulti("idac/Alkyne-Water.csv", modelClass));
//		experiments.add(new IDACExperimentsMulti("idac/AlkylHalide-Water.csv", modelClass));
//		experiments.add(new IDACExperimentsMulti("idac/Aromatic-Water.csv", modelClass));
//		experiments.add(new IDACExperimentsMulti("idac/ArylHalide-Water.csv", modelClass));
//		experiments.add(new IDACExperimentsMulti("idac/MultiringAromatics-Water.csv", modelClass));
////		experiments.add(new IDACExperimentsMulti("idac/CarboxilicAcid-Water.csv", modelClass)); //
//		experiments.add(new IDACExperimentsMulti("idac/CycloAlkane-Water.csv", modelClass));
////		experiments.add(new IDACExperimentsMulti("idac/CycloAlkene-Water.csv", modelClass));
////		experiments.add(new IDACExperimentsMulti("idac/Ether-Water.csv", modelClass)); //
////		experiments.add(new IDACExperimentsMulti("idac/Ester-Water.csv", modelClass));
//		experiments.add(new IDACExperimentsMulti("idac/Ketone-Water.csv", modelClass));
//		experiments.add(new IDACExperimentsMulti("idac/Water.csv", modelClass));
		
//		experiments.add(new IDACExperimentsMulti("idac/Alcohol-Alkane.csv", modelClass));
//		experiments.add(new IDACExperimentsMulti("idac/Alcohol-CycloAlkane.csv", modelClass));
//
//		experiments.add(new IDACExperimentsMulti("idac/Alkane-Alcohol.csv", modelClass));
////		experiments.add(new IDACExperimentsMulti("idac/Alkane-Alkane.csv", modelClass));
////		experiments.add(new IDACExperimentsMulti("idac/Alkane-AlkylHalide.csv", modelClass));
////		experiments.add(new IDACExperimentsMulti("idac/Alkane-Amine.csv", modelClass)); //
//		experiments.add(new IDACExperimentsMulti("idac/Alkane-CarboxilicAcid.csv", modelClass)); //
////		experiments.add(new IDACExperimentsMulti("idac/Alkane-Ketone.csv", modelClass)); //
//		experiments.add(new IDACExperimentsMulti("idac/Alkane-Phenol.csv", modelClass));
////		
////		experiments.add(new IDACExperimentsMulti("idac/Alkene-Amine.csv", modelClass)); //
//////		
////		experiments.add(new IDACExperimentsMulti("idac/AlkylHalide-Alkane.csv", modelClass));
////		experiments.add(new IDACExperimentsMulti("idac/Amine-Alkane.csv", modelClass)); //
////		experiments.add(new IDACExperimentsMulti("idac/Aromatic-Alkane.csv", modelClass));
////
//		experiments.add(new IDACExperimentsMulti("idac/CycloAlkane-Alcohol.csv", modelClass));
////		experiments.add(new IDACExperimentsMulti("idac/CycloAlkane-AlkylHalide.csv", modelClass));
////		experiments.add(new IDACExperimentsMulti("idac/CycloAlkane-Amine.csv", modelClass)); //
//		experiments.add(new IDACExperimentsMulti("idac/CycloAlkane-CarboxilicAcid.csv", modelClass)); //
//		experiments.add(new IDACExperimentsMulti("idac/CycloAlkane-Phenol.csv", modelClass));
//
//		//		
//		experiments.add(new IDACExperimentsMulti("idac/CarboxilicAcid-Alkane.csv", modelClass)); //
////		experiments.add(new IDACExperimentsMulti("idac/CarboxilicAcid-CycloAlkane.csv", modelClass)); //
//		
//
////		experiments.add(new IDACExperimentsMulti("idac/CycloAlkene-Amine.csv", modelClass)); //
//
//		experiments.add(new IDACExperimentsMulti("idac/Ketone-Alcohol.csv", modelClass)); //
//		experiments.add(new IDACExperimentsMulti("idac/Ketone-Alkane.csv", modelClass)); //
		
		// or just the two main groups
		experiments.add(new IDACExperimentsMulti("idac/nonHB.csv", modelClass));
//		experiments.add(new IDACExperimentsMulti("idac/nonaqueous.csv", modelClass));
//		experiments.add(new IDACExperimentsMulti("idac/glycerol.csv", modelClass));
//		experiments.add(new IDACExperimentsMulti("idac/aqueous.csv", modelClass));
	}

	public boolean getBounds(double[] xl, double[] xu) {
		getCurrent(xl);
		for (int i = 0; i < xu.length; i++) {
			double x = xl[i];
			xl[i] = x - 0.98*x;
			xu[i] = x + 0.98*x;
		}
		return true;
	}

	public int getNumberOfVariables() {
		return getNumberOfPars();
	}
	public double objectiveFunction(double[] x) {
		try {
			return cost(x);
		} catch (CostException e) {
			return 1000; // return a high cost to avoid this point
		}
	}
	public int getNumberOfPars(){
		return (ndesc*ndesc+ndesc)/2 + 1;
//		return 3;
	}
	public void getCurrent(double [] pars){
		int i=0;
		COSMOSACMulti cosmo = (COSMOSACMulti) experiments.get(0).getModels().get(0);

		pars[i++] = cosmo.getBeta(0);
//		pars[i++] = cosmo.getAnorm();
//		pars[i++] = cosmo.getFpol();
		for (int k = 0; k < ndesc; k++) {
			for (int l = k; l < ndesc; l++) {
				pars[i++] = cosmo.getFpol(k,l);
			}
		}
	}

	public double cost(double[] pars) throws CostException {
		double cost = 0;
		int NP = 0;

		for(IDACExperimentsMulti exp : experiments){
//			exp.setPrintGammas(true);
			for (int j = 0; j < exp.getModels().size(); j++) {
				if(!exp.getValid().get(j))
					continue;
				
				COSMOSACMulti cosmo = exp.getModels().get(j);
				
				int i=0;
				cosmo.setBeta(pars[i++]);
//				cosmo.setAnorm(pars[i++]);
//				cosmo.setFpol(pars[i++]);
				for (int k = 0; k < ndesc; k++) {
					for (int l = k; l < ndesc; l++) {
						cosmo.setFpol(k,l, pars[i++]);
					}
				}
				
				// update some internal variables
				cosmo.parametersChanged();
			}

			try {
				exp.calcDeviations();
			} catch (Exception e) {
				throw new CostException(e);
			}
			// avoid NaN
			if(exp.getNP()>0 && exp.getAARD()==exp.getAARD()){
				NP += exp.getNP();
				cost += exp.getAARD() * exp.getNP();
			}
		}
		cost /= NP;

		System.out.println("PARS: ");
		for (int i = 0; i < pars.length; i++) {
			System.out.print(pars[i] + " ");
		}
//		System.out.println(" COST:" + cost);
		System.out.println(" COST:" + cost + " NP:" + NP);
		return cost;
	}

	class Checker implements ConvergenceChecker {
		public boolean converged(PointCostPair[] pair) {
			double sum = 0;
			for (int i = 0; i < pair.length; i++) {
				sum += Math.abs(pair[0].getCost() - pair[i].getCost());
			}
			return sum/pair[0].getCost() < 1e-4;
		}
	}

	public static void main(String[] args) {
		IDAC_EstMulti_Fpol est;
		try {
			est = new IDAC_EstMulti_Fpol();
		} catch (Exception e1) {
			e1.printStackTrace();
			return;
		}
		
		double [] x0 = new double[est.getNumberOfPars()];
		double [] x1 = new double[est.getNumberOfPars()];

		est.getCurrent(x0);

		DirectSearchOptimizer optimizer;		
		optimizer = new NelderMead();
		//		optimizer = new MultiDirectional();
		for (int i = 0; i < x1.length; i++) {
			x1[i] = 0.98 * x0[i];
		}
		try {
			optimizer.minimize(est, 2000, est.new Checker(), x0, x1);
			x1 = optimizer.getMinima()[0].getPoint();
		} catch (Exception e) {
			e.printStackTrace();
		}
		System.out.println("Solution found: ");
		for (int i = 0; i < x1.length; i++) {
			System.out.print(x1[i] + " ");
		}
		System.out.println(" COST:" + optimizer.getMinima()[0].getCost());

//		DiRect direct = new DiRect(est);
//		direct.setMinSize(1e-6);
////		direct.setTransform(new PaloschiTransformation());
//		direct.setMaxEvals(2000);
//		direct.optimize();
//		x1 = direct.getMinValue();
//		System.out.println("Solution found: ");
//		for (int i = 0; i < x1.length; i++) {
//			System.out.print(x1[i] + " ");
//		}
//		est.objectiveFunction(direct.getMinValue());
//		System.out.println(" COST:" + direct.getMinObjective());
		
		System.out.println("\nStarted from: ");
		for (int i = 0; i < x0.length; i++) {
			System.out.print(x0[i] + " ");
		}
		System.out.println();
		
		// print a model friend version of the solution
		COSMOSACMulti cosmo = est.experiments.get(0).getModels().get(0);
		System.out.println("setBeta(" + cosmo.getBeta() + ");");
		System.out.println("setBeta(1, " + cosmo.getBeta(1) + ");");
		System.out.println("setBeta(2, " + cosmo.getBeta(2) + ");");
		System.out.println("setCHB(" + cosmo.getCHB() + ");");
		System.out.println("setCHB(0, 0, " + cosmo.getCHB(0,0) + ");");
		System.out.println("setCHB(0, 1, " + cosmo.getCHB(0,1) + ");");
		System.out.println("setCHB(0, 2, " + cosmo.getCHB(0,2) + ");");
		System.out.println("setCHB(1, 1, " + cosmo.getCHB(1,1) + ");");
		System.out.println("setCHB(1, 2, " + cosmo.getCHB(1,2) + ");");
		System.out.println("setCHB(2, 2, " + cosmo.getCHB(2,2) + ");");
		System.out.println("setSigmaHB(" + cosmo.getSigmaHB() + ");");
		System.out.println("setSigmaHB2(" + cosmo.getSigmaHB2() + ");");
//		System.out.println("setSigmaHB3(" + cosmo.getSigmaHB3() + ");");
		for (int i = 0; i < ndesc; i++) {
			for (int j = i; j < ndesc; j++) {
				System.out.println("setFpol(" + i + ", " + j + ", " + cosmo.getFpol(i,j) + ");");
			}
		}
		System.out.println("setIgnoreSG(" + cosmo.isIgnoreSG() + ");");
		System.out.println("setCoord(" + cosmo.getCoord() + ");");
		System.out.println("setAnorm(" + cosmo.getAnorm() + ");");
		System.out.println("setVnorm(" + cosmo.getVnorm() + ");");
	}
}