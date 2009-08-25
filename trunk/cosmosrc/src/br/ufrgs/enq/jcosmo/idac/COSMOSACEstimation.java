package br.ufrgs.enq.jcosmo.idac;

import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math.estimation.EstimatedParameter;
import org.apache.commons.math.estimation.EstimationException;
import org.apache.commons.math.estimation.Estimator;
import org.apache.commons.math.estimation.LevenbergMarquardtEstimator;
import org.apache.commons.math.estimation.SimpleEstimationProblem;
import org.apache.commons.math.estimation.WeightedMeasurement;

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
		// add the estimated parameters    ("Name",      default value,      fixed (do not estimate) )
		addParameter(new EstimatedParameter("AEffPrime", COSMOSAC.AEFFPRIME, false));
		addParameter(new EstimatedParameter("Coord", COSMOSAC.COORD, false));
		addParameter(new EstimatedParameter("Vnorm", COSMOSAC.VNORM, false));
		addParameter(new EstimatedParameter("Anorm", COSMOSAC.ANORM, false));
		addParameter(new EstimatedParameter("CHB", COSMOSAC.CHB, false));
		addParameter(new EstimatedParameter("SigmaHB", COSMOSAC.SIGMAHB, false));
		addParameter(new EstimatedParameter("Epsilon", COSMOSAC.EPSILON, false));
		
		
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
		
		addAllMeasurements();
	}
	
	void addAllMeasurements(){
		for(IDACExperiments e : experiments){
			List<Double> measures = e.getMeasures();
			List<COSMOSAC> models = e.getModels();
			for (int i = 0; i < measures.size(); i++) {
				addMeasurement(new Measure(measures.get(i), models.get(i)));
			}
		}
	}

	public static void main(String[] args) throws EstimationException {
		COSMOSACEstimation est;
		try {
			est = new COSMOSACEstimation();
		} catch (Exception e1) {
			e1.printStackTrace();
			return;
		}

		Estimator solver = new LevenbergMarquardtEstimator();
//		Estimator solver = new GaussNewtonEstimator();
		
		solver.estimate(est);
		double[] errors = solver.guessParametersErrors(est);
		
		EstimatedParameter[] pars = est.getUnboundParameters();
		
		for (int j = 0; j < errors.length; j++) {
			EstimatedParameter p = pars[j];
			System.out.println(p.getName() + "= " + p.getEstimate() + " +- " + errors[j]);
		}
	}
}
