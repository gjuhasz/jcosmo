package br.ufrgs.enq.jcosmo.idac;

import java.awt.BasicStroke;
import java.awt.BorderLayout;
import java.awt.Color;
import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.Vector;

import javax.swing.JFrame;
import javax.swing.JPanel;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.labels.XYToolTipGenerator;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import org.jfree.data.Range;
import org.jfree.data.xy.XYDataset;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

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
public class IDACDiagonal extends JFrame implements XYToolTipGenerator{
	private static final long serialVersionUID = 1L;

	private String modelClass;
	private int NP=0;
	private int start;
	private Boolean useAll=true;
	
	List<COSMOSAC> models = new ArrayList<COSMOSAC>();
	
	List<Double> gammaInfs = new ArrayList<Double>();
	List<Double> temperatures = new ArrayList<Double>();
	
	List<XYSeries> points = new ArrayList<XYSeries>();
	List<List<String>> tps = new Vector<List<String>>();

	List<Color> color = new ArrayList<Color>();
	
	public void addIDACExperiments(String filename, String modelClass, Color cor) throws Exception {
		this.modelClass = modelClass;
		color.add(cor);
		
		start = models.size();
		loadContents(filename);
		addPoints(filename);
	}
	
	public void addIDACExperiments(String filename, String modelClass) throws Exception {
		addIDACExperiments(filename,modelClass, null);
	}
	
	public void loadContents(String filename) throws Exception{
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
			if (modelClass != COSMOSAC.class.getName() & useAll==true){
				comps[0] = db.getComp("water");
				comps[0].name = compNames[0].toUpperCase();
				comps[1] = db.getComp("water");
				comps[1].name = compNames[1].toUpperCase();
			} 
			else {
				comps[0] = db.getComp(compNames[0].replace(' ', '-'));
				comps[1] = db.getComp(compNames[1].replace(' ', '-'));
				if(comps[0] == null){
					System.err.println("Component " + compNames[0] + " not found, it will be ignored.");
					continue;
				}
				if(comps[1] == null){
					System.err.println("Component " + compNames[1] + " not found, it will be ignored.");
					continue;
				}
			}
			
			try{
				model.setComponents(comps);
			}
			catch (Exception e) {
				System.err.println(e.getMessage() + ", it will be ignored.");
				continue;
			}
			
			models.add(model);
			temperatures.add(T);
			gammaInfs.add(gammaInf);
			NP++;
		}
		reader.close();
	}
	
	public void addPoints(String filename) throws Exception{
		
		double z[] = new double[2];
		double lnGamma[] = new double[2];
		double T;
		// infinite dilution activity calculation
		z[0] = 0; z[1] = 1-z[0];
		File file = new File(filename);
		XYSeries point = new XYSeries(file.getName().substring(0, file.getName().lastIndexOf('.')), false);
		List<String> tpNames = new Vector<String>();
		
		for (int i = start; i < models.size(); i++) {
			COSMOSAC model = models.get(i);
			T = temperatures.get(i);
			double gammaInf = gammaInfs.get(i);
			
			model.setTemperature(T);
			model.setComposition(z);
			model.activityCoefficientLn(lnGamma, 0);
						
			point.add(Math.log(gammaInf), lnGamma[0]);
			tpNames.add(model.getComps()[0].name + "/" + model.getComps()[1].name +
					", T=" + T + ", " + gammaInf +  " (" + Math.log(gammaInf) + ", " + lnGamma[0] + ")");
		}
		points.add(point);
		tps.add(tpNames);
	}
	
	public void showPlot(){
	
		XYSeriesCollection dataset = new XYSeriesCollection();

		XYPlot plot;
		JFreeChart chart = ChartFactory.createXYLineChart(null, 
				"Logarithm of Experimental IDAC", "Logarithm of Model IDAC", dataset,
				PlotOrientation.VERTICAL, true, true, false);
		plot = (XYPlot) chart.getPlot();		
		plot.getDomainAxis().setRange(new Range(-2, 26));
		plot.getRangeAxis().setRange(new Range(-2, 26));
		
		plot.setBackgroundPaint(Color.black);
	
		XYLineAndShapeRenderer r = (XYLineAndShapeRenderer) plot.getRenderer();
		int i;
		for(i=0; i<points.size(); i++){
			dataset.addSeries(points.get(i));
			r.setSeriesLinesVisible(i, false);
			r.setSeriesShapesVisible(i, true);
			r.setSeriesToolTipGenerator(i, this);
		}
		
		BasicStroke dashed = new BasicStroke(1.0f, BasicStroke.CAP_ROUND, BasicStroke.JOIN_ROUND,  
				1.0f, new float[] { 10.0f, 6.0f }, 0.0f);
		
		XYSeries diag = new XYSeries(NP);
		diag.add(-2, -2);
		diag.add(26, 26);
		dataset.addSeries(diag);
		r.setSeriesStroke(i, dashed);
		r.setSeriesLinesVisible(i++, true);
		
		XYSeries erroPos = new XYSeries("-1 ln unit");
		erroPos.add(-2, -1);
		erroPos.add(26, 27);
		dataset.addSeries(erroPos);
		r.setSeriesStroke(i, dashed);
		r.setSeriesLinesVisible(i++, true);
		
		XYSeries erroNeg = new XYSeries("+1 ln unit");
		erroNeg.add(-2, -3);
		erroNeg.add(26, 25);
		dataset.addSeries(erroNeg);
		r.setSeriesStroke(i, dashed);
		r.setSeriesLinesVisible(i++, true);
		
		ChartPanel chartPanel = new ChartPanel(chart);
		JPanel jp = new JPanel(new BorderLayout());
		jp.add(chartPanel);
		
		add(jp, BorderLayout.CENTER);
		setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		setSize(400,500);
		setVisible(true);	
	}
	
	public void setUseAll(Boolean useAll) {
		this.useAll = useAll;
	}

	public String generateToolTip(XYDataset dataset, int series, int item) {
		return tps.get(series).get(item);
	}
}
