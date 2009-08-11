package br.ufrgs.enq.jcosmo.idac;

import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math.estimation.EstimatedParameter;
import org.apache.commons.math.estimation.EstimationException;
import org.apache.commons.math.estimation.SimpleEstimationProblem;
import org.apache.commons.math.estimation.WeightedMeasurement;
import org.apache.commons.math.optimization.ConvergenceChecker;
import org.apache.commons.math.optimization.CostException;
import org.apache.commons.math.optimization.CostFunction;
import org.apache.commons.math.optimization.DirectSearchOptimizer;
import org.apache.commons.math.optimization.NelderMead;
import org.apache.commons.math.optimization.PointCostPair;

import br.ufrgs.enq.direct.DiRect;
import br.ufrgs.enq.direct.ObjectiveFunction;
import br.ufrgs.enq.direct.PaloschiTransformation;
import br.ufrgs.enq.jcosmo.COSMOPAC;
import br.ufrgs.enq.jcosmo.COSMOSAC;

public class COSMOSACEstimation extends SimpleEstimationProblem {
	
	List<IDACExperiments> experiments;
	
	class Measure extends WeightedMeasurement {
		COSMOSAC model;
		double lnGamma[] = new double[2];
		
		public Measure(double value, COSMOSAC model) {
			super(1.0, value);
			this.model = model;
		}

		@Override
		public double getPartial(EstimatedParameter p) {
			double oldy = getTheoreticalValue();
			double oldPar = p.getEstimate();
			double delta = 1e-4;
			
			p.setEstimate(oldPar + delta);
			double newy = getTheoreticalValue();
			
			p.setEstimate(oldPar);
			
			return (newy - oldy)/delta;
		}
		@Override
		public double getTheoreticalValue() {
			for(EstimatedParameter p : getUnboundParameters()){
				if(p.getName().equals("AEffPrime"))
					model.setAEffPrime(p.getEstimate());
				else if(p.getName().equals("Coord"))
					model.setCoord(p.getEstimate());
				else if(p.getName().equals("Vnorm"))
					model.setVnorm(p.getEstimate());
				else if(p.getName().equals("Anorm"))
					model.setAnorm(p.getEstimate());
				else if(p.getName().equals("CHB"))
					model.setCHB(p.getEstimate());
				else if(p.getName().equals("SigmaHB"))
					model.setSigmaHB(p.getEstimate());
				else if(p.getName().equals("Epsilon"))
					model.setEpsilon(p.getEstimate());
			}
			model.parametersChanged();
			
			model.activityCoefficient(lnGamma);
			return lnGamma[0];
		}
		
	}

	public COSMOSACEstimation() throws Exception {
//		String modelClass = COSMOSAC.class.getName();
//		String modelClass = COSMOSACGAMESS.class.getName();
		String modelClass = COSMOPAC.class.getName();

		experiments = new ArrayList<IDACExperiments>();

//		File dir = new File("idac");
//		String[] files = dir.list();
//		for(String f : files){
//			if(f.endsWith(".csv") 
//					&& !f.contains("Amine")
//					&& !f.contains("Acid")
//					&& !f.contains("Heavy")
//					){
//				IDACExperiments exp = new IDACExperiments(dir.getName() + '/' + f, modelClass, false);
//				experiments.add(exp);
//			}
//		}

		experiments.add(new IDACExperiments("idac/Alcohol-Water.csv", modelClass));
		experiments.add(new IDACExperiments("idac/Aldehyde-Water.csv", modelClass));
		experiments.add(new IDACExperiments("idac/Alkane-Water.csv", modelClass));
		experiments.add(new IDACExperiments("idac/AlkylHalide-Water.csv", modelClass));
		experiments.add(new IDACExperiments("idac/Aromatic-Water.csv", modelClass));
		experiments.add(new IDACExperiments("idac/MultiringAromatics-Water.csv", modelClass));
//		experiments.add(new IDACExperiments("idac/CarboxilicAcid-Water.csv", modelClass));//
		experiments.add(new IDACExperiments("idac/CycloAlkane-Water.csv", modelClass));
		experiments.add(new IDACExperiments("idac/Ketone-Water.csv", modelClass));
		experiments.add(new IDACExperiments("idac/Water.csv", modelClass));
		
//		experiments.add(new IDACExperiments("idac/aqueous.csv", modelClass));
				
		experiments.add(new IDACExperiments("idac/Alcohol-Alkane.csv", modelClass));
		experiments.add(new IDACExperiments("idac/Alcohol-CycloAlkane.csv", modelClass));

		experiments.add(new IDACExperiments("idac/Alkane-Alcohol.csv", modelClass));
		experiments.add(new IDACExperiments("idac/Alkane-AlkylHalide.csv", modelClass));
		experiments.add(new IDACExperiments("idac/Alkane-Amine.csv", modelClass));
//		experiments.add(new IDACExperiments("idac/Alkane-CarboxilicAcid.csv", modelClass));//
		experiments.add(new IDACExperiments("idac/Alkane-Ketone.csv", modelClass));
		experiments.add(new IDACExperiments("idac/Alkane-Phenol.csv", modelClass));
		
		experiments.add(new IDACExperiments("idac/Alkene-Amine.csv", modelClass));
		
		experiments.add(new IDACExperiments("idac/AlkylHalide-Alkane.csv", modelClass));
//		experiments.add(new IDACExperiments("idac/Amine-Alkane.csv", modelClass));//
		experiments.add(new IDACExperiments("idac/Aromatic-Alkane.csv", modelClass));
		
		experiments.add(new IDACExperiments("idac/CycloAlkane-Alcohol.csv", modelClass));
		experiments.add(new IDACExperiments("idac/CycloAlkane-AlkylHalide.csv", modelClass));
		experiments.add(new IDACExperiments("idac/CycloAlkane-Amine.csv", modelClass));
//		experiments.add(new IDACExperiments("idac/CycloAlkane-CarboxilicAcid.csv", modelClass));//
		experiments.add(new IDACExperiments("idac/CycloAlkane-Phenol.csv", modelClass));
		
		experiments.add(new IDACExperiments("idac/CycloAlkene-Amine.csv", modelClass));

//		experiments.add(new IDACExperiments("idac/Ketone-Alcohol.csv", modelClass));//
		experiments.add(new IDACExperiments("idac/Ketone-Alkane.csv", modelClass));
		
//		experiments.add(new IDACExperiments("idac/nonaqueous.csv", modelClass));
//		experiments.add(new IDACExperiments("idac/nonaqueous2.csv", modelClass));
		
		experiments.add(new IDACExperiments("idac/Alkane-Alkane.csv", modelClass));
	}

	public boolean getBounds(double[] xl, double[] xu) {
		getCurrent(xl);
		for (int i = 0; i < xu.length; i++) {
			double x = xl[i];
			xl[i] = x - 0.8*x;
			xu[i] = x + 0.8*x;
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
		return 4;
	}
	public void getCurrent(double [] pars){
		int i=0;
		COSMOSAC cosmo = (COSMOSAC) experiments.get(0).getModels().get(0);

//		pars[i++] = cosmo.getAEffPrime();
//		pars[i++] = cosmo.getCoord();
//		pars[i++] = cosmo.getVnorm();
		pars[i++] = cosmo.getAnorm();
		pars[i++] = cosmo.getCHB();
		pars[i++] = cosmo.getSigmaHB();
		pars[i++] = cosmo.getEpsilon();
	}

	public double cost(double[] pars) throws CostException {
		double cost = 0;
		int NP = 0;

		for(IDACExperiments exp : experiments){
			for(COSMOSAC model : exp.getModels()){

				COSMOSAC cosmo = (COSMOSAC) model;

				int i=0;
//				cosmo.setAEffPrime(pars[i++]);
//				cosmo.setCoord(pars[i++]);
//				cosmo.setVnorm(pars[i++]);
				cosmo.setAnorm(pars[i++]);
				cosmo.setCHB(pars[i++]);
				cosmo.setSigmaHB(pars[i++]);
				cosmo.setEpsilon(pars[i++]);
				
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
			return Math.abs(pair[0].getCost() - pair[1].getCost())/pair[0].getCost() < 1e-3;
		}
	}

	public static void main(String[] args) {
		COSMOSACEstimation est;
		try {
			est = new COSMOSACEstimation();
		} catch (Exception e1) {
			e1.printStackTrace();
			return;
		}

		double [] x0 = new double[est.getNumberOfPars()];
		double [] x1 = new double[est.getNumberOfPars()];

		est.getCurrent(x0);

//		DirectSearchOptimizer optimizer;		
//		optimizer = new NelderMead();
//		//		optimizer = new MultiDirectional();
//		for (int i = 0; i < x1.length; i++) {
//			x1[i] = 0.8 * x0[i];
//		}
//		try {
//			optimizer.minimize(est, 2000, est.new Checker(), x0, x1);
//			x1 = optimizer.getMinima()[0].getPoint();
//		} catch (Exception e) {
//			e.printStackTrace();
//		}
//		System.out.println("Solution found: ");
//		for (int i = 0; i < x1.length; i++) {
//			System.out.print(x1[i] + " ");
//		}
//		System.out.println(" COST:" + optimizer.getMinima()[0].getCost());

		System.out.println("\nStarted from: ");
		for (int i = 0; i < x0.length; i++) {
			System.out.print(x0[i] + " ");
		}
	}
}
