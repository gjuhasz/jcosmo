package br.ufrgs.enq.jcosmo.test;

import java.awt.BorderLayout;
import java.awt.Color;
import java.util.ArrayList;
import java.util.List;

import javax.swing.JFrame;
import javax.swing.JPanel;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYSplineRenderer;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

import br.ufrgs.enq.jcosmo.COSMOPAC;
import br.ufrgs.enq.jcosmo.COSMOSAC;
import br.ufrgs.enq.jcosmo.COSMOSACCompound;
import br.ufrgs.enq.jcosmo.COSMOSACDataBase;

import com.csvreader.CsvReader;

/**
 * Class representing an set of infinite dilution activity coefficient (IDAC) experiments.
 * 
 * @author rafael e renan
 *
 */
public class Diagonal extends JFrame {
	private static final long serialVersionUID = 1L;

	String filename;

	private String modelClass;
	private int NP=0;
	
	List<COSMOSAC> models = new ArrayList<COSMOSAC>();
	
	List<Double> gammaInfs = new ArrayList<Double>();
	List<Double> temperatures = new ArrayList<Double>();
	
	List<XYSeries> points = new ArrayList<XYSeries>();
	List<Color> color = new ArrayList<Color>();
//	static XYSeries point = new XYSeries(" ");
	
	public void addIDACExperiments(String filename, String modelClass, Color cor) throws Exception {
		this.filename = filename;
		this.modelClass = modelClass;
		color.add(cor);
		
		loadContents();
		addPoint();
	}
	
	public void addIDACExperiments(String filename, String modelClass) throws Exception {
		addIDACExperiments(filename,modelClass, null);
	}
	
	public void loadContents() throws Exception{
		COSMOSAC model = null;
		
		CsvReader reader = new CsvReader(filename);
		
		// skip the headers
		reader.readHeaders();
		
		COSMOSACDataBase db = COSMOSACDataBase.getInstance();
		
		while(reader.readRecord()){
			model = (COSMOSAC) Class.forName(modelClass).newInstance();
			
			String mixtureString = reader.get(0);
			double T, gammaInf;
			try{
				T = Double.parseDouble(reader.get(1));
				gammaInf = Double.parseDouble(reader.get(2));
			}
			catch (NumberFormatException e) {
				System.err.println(e.toString() + ", line " + (NP + 1));
				continue;
			}
			
			String[] compNames = mixtureString.split("/");
			if(compNames.length != 2) {
				System.err.println("Typo on mixture(" + mixtureString + ", line " + (NP + 1));
				continue;
			}
			
			COSMOSACCompound comps[] = new COSMOSACCompound[2];
			if (modelClass == COSMOPAC.class.getName()){
				comps[0] = db.getComp("water");
				comps[0].name = compNames[0].toUpperCase();
				comps[1] = db.getComp("water");
				comps[1].name = compNames[1].toUpperCase();
			} 
			else {
				comps[0] = db.getComp(compNames[0].replace(' ', '-'));
				comps[1] = db.getComp(compNames[1].replace(' ', '-'));
				if(comps[0] == null)
					throw new IllegalArgumentException("Component " + compNames[0] + " not found.");
				if(comps[1] == null)
					throw new IllegalArgumentException("Component " + compNames[1] + " not found.");
			}
			
			model.setComponents(comps);
			
			models.add(model);
			temperatures.add(T);
			gammaInfs.add(gammaInf);
		}
		reader.close();
	}
	
	public void addPoint() throws Exception{
		
		double z[] = new double[2];
		double lnGamma[] = new double[2];
		double T;
		// infinite dilution activity calculation
		z[0] = 0; z[1] = 1-z[0];
		XYSeries point = new XYSeries(filename);
		
		for (int i = 0; i < models.size(); i++) {
			COSMOSAC model = models.get(i);
			T = temperatures.get(i);
			double gammaInf = gammaInfs.get(i);
			
			model.setTemperature(T);
			model.setComposition(z);
			model.activityCoefficientLn(lnGamma, 0);
						
			point.add(lnGamma[0], Math.log(gammaInf));
			NP++;
		}
		points.add(point);
	}
	
	public void PlotDig(){
				
		XYPlot plot;
		JFreeChart chart = ChartFactory.createXYLineChart(null, 
				"calc", "exp", null, PlotOrientation.VERTICAL, true, true, false);
		plot = (XYPlot) chart.getPlot();		
//		plot.getDomainAxis().setRange(new Range(-1, 24));
//		plot.getRangeAxis().setRange(new Range(-1, 24));
	
		XYSeriesCollection[] dataset = new XYSeriesCollection [points.size()];
		int i;
		for(i=0; i<points.size(); i++){
			XYSplineRenderer r = new XYSplineRenderer();
			r.setBaseLinesVisible(false);
			dataset[i] = new XYSeriesCollection();
			dataset[i].addSeries(points.get(i));
			plot.setDataset(i, dataset[i]);
//			r.setBaseFillPaint(color.get(i));
			plot.setRenderer(i, r);
		}		
		
//		for(i=0; i<points.size(); i++){
//			plot.setDataset(i, dataset[i]);
//		}

//		XYSeriesCollection dataset = new XYSeriesCollection();
//		dataset.addSeries(points.get(1));
//		dataset.addSeries(points.get(2));
//		plot.setDataset(0, dataset);
//		XYSplineRenderer r = new XYSplineRenderer();
//		r.setBaseLinesVisible(false);
//		plot.setRenderer(0, r);
		
		XYSeriesCollection dataset2 = new XYSeriesCollection();
		XYSeries diag = new XYSeries("diagonal");
		diag.add(-1, -1);
		diag.add(24, 24);
		dataset2.addSeries(diag);
		
		plot.setDataset(1, dataset2);
		XYSplineRenderer r2 = new XYSplineRenderer();
		r2.setBaseShapesVisible(false);
		plot.setRenderer(1, r2);
				
		ChartPanel chartPanel = new ChartPanel(chart);
		JPanel jp = new JPanel(new BorderLayout());
		jp.add(chartPanel);
		
		add(jp, BorderLayout.CENTER);
		setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		setSize(400,500);
		setVisible(true);	
	}
		
	public int getNP() {
		return NP;
	}
	
	public List<COSMOSAC> getModels() {
		return models;
	}
	
	
	public static void main (String[] args) throws Exception{
		String modelClass = COSMOPAC.class.getName();
		
		Diagonal dig = new Diagonal();
		dig.addIDACExperiments("idac/Aldehyde-Water.csv", modelClass);
		dig.addIDACExperiments("idac/Alkane-Water.csv", modelClass);
		dig.addIDACExperiments("idac/Aromatic-Water.csv", modelClass);
		dig.addIDACExperiments("idac/ChloroAlkane-Water.csv", modelClass);
		dig.addIDACExperiments("idac/CycloAlkane-Water.csv", modelClass);
		dig.addIDACExperiments("idac/Water-Alkane.csv", modelClass);
		dig.addIDACExperiments("idac/Water-CycloAlkane.csv", modelClass);
		dig.addIDACExperiments("idac/Water-Alcohol.csv", modelClass);
		dig.addIDACExperiments("idac/Water-HeavyAlcohol.csv", modelClass);
		
//		dig.addIDACExperiments("idac/Alkane-Alcohol.csv", modelClass);
//		dig.addIDACExperiments("idac/Alcohol-Alkane.csv", modelClass);
//
//		dig.addIDACExperiments("idac/Alcohol-HeavyAlkane.csv", modelClass);
//		dig.addIDACExperiments("idac/HeavyAlkane-Alcohol.csv", modelClass);
//
//		dig.addIDACExperiments("idac/Alcohol-CycloAlkane.csv", modelClass);
//		dig.addIDACExperiments("idac/CycloAlkane-Alcohol.csv", modelClass);
//
//		dig.addIDACExperiments("idac/Aromatic-Alkane.csv", modelClass);
//		dig.addIDACExperiments("idac/CycloAlkane-Phenol.csv", modelClass);
//		dig.addIDACExperiments("idac/Alkane-Phenol.csv", modelClass);
//
//		dig.addIDACExperiments("idac/Alkane-AlkylHalide.csv", modelClass);
//		dig.addIDACExperiments("idac/AlkylHalide-Alkane.csv", modelClass);
		
//		dig.addIDACExperiments("idac/Ketone-Alkane.csv", modelClass, Color.RED);
//		dig.addIDACExperiments("idac/Ketone-Alcohol.csv", modelClass);
//		dig.addIDACExperiments("idac/CarboxilicAcid-Alkane.csv", modelClass);
//		dig.addIDACExperiments("idac/Alkane-CarboxilicAcid.csv", modelClass);
//		dig.addIDACExperiments("idac/Ketone-Water.csv", modelClass);
		
		dig.PlotDig();	
	}
}