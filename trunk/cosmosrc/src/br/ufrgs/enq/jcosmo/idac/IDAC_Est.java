package br.ufrgs.enq.jcosmo.idac;

import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math.optimization.ConvergenceChecker;
import org.apache.commons.math.optimization.CostException;
import org.apache.commons.math.optimization.CostFunction;
import org.apache.commons.math.optimization.DirectSearchOptimizer;
import org.apache.commons.math.optimization.NelderMead;
import org.apache.commons.math.optimization.PointCostPair;

import br.ufrgs.enq.direct.DiRect;
import br.ufrgs.enq.direct.ObjectiveFunction;
import br.ufrgs.enq.jcosmo.COSMOPAC;
import br.ufrgs.enq.jcosmo.COSMOSAC;
import br.ufrgs.enq.jcosmo.COSMOSAC_G;
import br.ufrgs.enq.jcosmo.PCMSAC;

public class IDAC_Est implements CostFunction, ObjectiveFunction {

	List<IDACExperiments> experiments;

	public IDAC_Est() throws Exception {
//		String modelClass = COSMOSAC.class.getName();
		String modelClass = COSMOSAC_G.class.getName();
//		String modelClass = COSMOPAC.class.getName();
//		String modelClass = PCMSAC.class.getName();

		experiments = new ArrayList<IDACExperiments>();

//		experiments.add(new IDACExperiments("idac/Alcohol-Water.csv", modelClass));
//		experiments.add(new IDACExperiments("idac/Aldehyde-Water.csv", modelClass));
//		experiments.add(new IDACExperiments("idac/Alkane-Water.csv", modelClass));
//		experiments.add(new IDACExperiments("idac/Alkene-Water.csv", modelClass));
//		experiments.add(new IDACExperiments("idac/Alkyne-Water.csv", modelClass));
//		experiments.add(new IDACExperiments("idac/AlkylHalide-Water.csv", modelClass));
//		experiments.add(new IDACExperiments("idac/Aromatic-Water.csv", modelClass));
//		experiments.add(new IDACExperiments("idac/ArylHalide-Water.csv", modelClass));
//		experiments.add(new IDACExperiments("idac/MultiringAromatics-Water.csv", modelClass));
////		experiments.add(new IDACExperiments("idac/CarboxilicAcid-Water.csv", modelClass)); //
//		experiments.add(new IDACExperiments("idac/CycloAlkane-Water.csv", modelClass));
//		experiments.add(new IDACExperiments("idac/CycloAlkene-Water.csv", modelClass));
//		experiments.add(new IDACExperiments("idac/Ether-Water.csv", modelClass)); //
//		experiments.add(new IDACExperiments("idac/Ester-Water.csv", modelClass));
//		experiments.add(new IDACExperiments("idac/Ketone-Water.csv", modelClass));
//		experiments.add(new IDACExperiments("idac/Water.csv", modelClass));

		//		
//		experiments.add(new IDACExperiments("idac/Alcohol-Alkane.csv", modelClass));
//		experiments.add(new IDACExperiments("idac/Alcohol-CycloAlkane.csv", modelClass));
//
//		experiments.add(new IDACExperiments("idac/Alkane-Alcohol.csv", modelClass));
//		experiments.add(new IDACExperiments("idac/Alkane-Alkane.csv", modelClass));
//		experiments.add(new IDACExperiments("idac/Alkane-AlkylHalide.csv", modelClass));
//		experiments.add(new IDACExperiments("idac/Alkane-Amine.csv", modelClass)); //
//		experiments.add(new IDACExperiments("idac/Alkane-CarboxilicAcid.csv", modelClass)); //
//		experiments.add(new IDACExperiments("idac/Alkane-Ketone.csv", modelClass)); //
//		experiments.add(new IDACExperiments("idac/Alkane-Phenol.csv", modelClass));
		
//		experiments.add(new IDACExperiments("idac/Alkene-Amine.csv", modelClass)); //
//		
//		experiments.add(new IDACExperiments("idac/AlkylHalide-Alkane.csv", modelClass));
//		experiments.add(new IDACExperiments("idac/Amine-Alkane.csv", modelClass)); //
//		experiments.add(new IDACExperiments("idac/Aromatic-Alkane.csv", modelClass));

//		experiments.add(new IDACExperiments("idac/CycloAlkane-Alcohol.csv", modelClass));
//		experiments.add(new IDACExperiments("idac/CycloAlkane-AlkylHalide.csv", modelClass));
//		experiments.add(new IDACExperiments("idac/CycloAlkane-Amine.csv", modelClass)); //
//		experiments.add(new IDACExperiments("idac/CycloAlkane-CarboxilicAcid.csv", modelClass)); //
//		experiments.add(new IDACExperiments("idac/CycloAlkane-Phenol.csv", modelClass));

		//		
////		experiments.add(new IDACExperiments("idac/CarboxilicAcid-Alkane.csv", modelClass)); //
////		experiments.add(new IDACExperiments("idac/CarboxilicAcid-CycloAlkane.csv", modelClass)); //
//		
//
////		experiments.add(new IDACExperiments("idac/CycloAlkene-Amine.csv", modelClass)); //
//
//		experiments.add(new IDACExperiments("idac/Ketone-Alcohol.csv", modelClass)); //
//		experiments.add(new IDACExperiments("idac/Ketone-Alkane.csv", modelClass)); //
		
		// or just the two main groups
		experiments.add(new IDACExperiments("idac/nonHB.csv", modelClass));
//		experiments.add(new IDACExperiments("idac/nonaqueous.csv", modelClass));
//		experiments.add(new IDACExperiments("idac/aqueous.csv", modelClass));
	}

	public boolean getBounds(double[] xl, double[] xu) {
		getCurrent(xl);
		for (int i = 0; i < xu.length; i++) {
			double x = xl[i];
			xl[i] = 0.1*x;
			xu[i] = 10*x;
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
//		return 6;
		return 2;
	}
	public void getCurrent(double [] pars){
		int i=0;
		COSMOSAC cosmo = (COSMOSAC) experiments.get(0).getModels().get(0);

		pars[i++] = cosmo.getBeta();
		pars[i++] = cosmo.getFpol();
//		pars[i++] = cosmo.getCHB();
//		pars[i++] = cosmo.getSigmaHB();
//		pars[i++] = cosmo.getSigmaHB2();
//		pars[i++] = cosmo.getSigmaHB3();
//		pars[i++] = cosmo.getCoord();
//		pars[i++] = cosmo.getAnorm();
//		pars[i++] = cosmo.getVnorm();
	}

	public double cost(double[] pars) throws CostException {
		double cost = 0;
		int NP = 0;

		for(IDACExperiments exp : experiments){
//			exp.setPrintGammas(true);
			for (int j = 0; j < exp.getModels().size(); j++) {
				if(!exp.getValid().get(j))
					continue;
				
				COSMOSAC cosmo = exp.getModels().get(j);
				
				int i=0;
				cosmo.setBeta(pars[i++]);
				cosmo.setFpol(pars[i++]);
//				cosmo.setCHB(pars[i++]);
//				cosmo.setSigmaHB(pars[i++]);
//				cosmo.setSigmaHB2(pars[i++]);
//				cosmo.setSigmaHB3(pars[i++]);
//				cosmo.setCoord(pars[i++]);
//				cosmo.setAnorm(pars[i++]);
//				cosmo.setVnorm(pars[i++]);
				
				// update some internal variables
				cosmo.parametersChanged();
			}

			try {
				exp.calcDeviations();
			} catch (Exception e) {
				throw new CostException(e);
			}
			if(exp.getNP()>0){
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
			return sum/pair[0].getCost() < 1e-5;
		}
	}

	public static void main(String[] args) {
		IDAC_Est est;
		try {
			est = new IDAC_Est();
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
			x1[i] = 0.8 * x0[i];
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
//		
//		System.out.println("\nStarted from: ");
//		for (int i = 0; i < x0.length; i++) {
//			System.out.print(x0[i] + " ");
//		}
//		System.out.println();
		
		// print a model friend version of the solution
		COSMOSAC cosmo = est.experiments.get(0).getModels().get(0);
		System.out.println("setBeta(" + cosmo.getBeta() + ");");
		System.out.println("setCHB(" + cosmo.getCHB() + ");");
		System.out.println("setSigmaHB(" + cosmo.getSigmaHB() + ");");
		System.out.println("setSigmaHB2(" + cosmo.getSigmaHB2() + ");");
		System.out.println("setSigmaHB3(" + cosmo.getSigmaHB3() + ");");
		System.out.println("setFpol(" + cosmo.getFpol() + ");");
		System.out.println("setIgnoreSG(" + cosmo.isIgnoreSG() + ");");
		System.out.println("setCoord(" + cosmo.getCoord() + ");");
		System.out.println("setAnorm(" + cosmo.getAnorm() + ");");
		System.out.println("setVnorm(" + cosmo.getVnorm() + ");");
	}
}
