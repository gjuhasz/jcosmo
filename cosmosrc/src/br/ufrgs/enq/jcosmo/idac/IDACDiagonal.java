package br.ufrgs.enq.jcosmo.idac;

import java.awt.BasicStroke;
import java.awt.BorderLayout;
import java.awt.Color;
import java.io.File;
import java.util.ArrayList;
import java.util.List;

import javax.swing.JFrame;
import javax.swing.JPanel;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import org.jfree.chart.renderer.xy.XYSplineRenderer;
import org.jfree.data.Range;
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
public class IDACDiagonal extends JFrame {
	private static final long serialVersionUID = 1L;

	private String modelClass;
	private int NP=0;
	private int start;
	private Boolean useAll=true;
	
	List<COSMOSAC> models = new ArrayList<COSMOSAC>();
	
	List<Double> gammaInfs = new ArrayList<Double>();
	List<Double> temperatures = new ArrayList<Double>();
	
	List<XYSeries> points = new ArrayList<XYSeries>();
	List<Color> color = new ArrayList<Color>();
//	static XYSeries point = new XYSeries(" ");
	
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
			if (modelClass == COSMOPAC.class.getName() & useAll==true){
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
		XYSeries point = new XYSeries(file.getName().substring(0, file.getName().lastIndexOf('.')));
		
		for (int i = start; i < models.size(); i++) {
			COSMOSAC model = models.get(i);
			T = temperatures.get(i);
			double gammaInf = gammaInfs.get(i);
			
			model.setTemperature(T);
			model.setComposition(z);
			model.activityCoefficientLn(lnGamma, 0);
						
			point.add(Math.log(gammaInf), lnGamma[0]);
		}
		points.add(point);
	}
	
	@SuppressWarnings("deprecation")
	public void showPlot(){
				
		XYPlot plot;
		JFreeChart chart = ChartFactory.createXYLineChart(null, 
				"Logarithm of Experimental IDAC", "Logarithm of Model IDAC", null, PlotOrientation.VERTICAL, true, true, false);
		plot = (XYPlot) chart.getPlot();		
		plot.getDomainAxis().setRange(new Range(-2, 26));
		plot.getRangeAxis().setRange(new Range(-2, 26));
		
		plot.setBackgroundPaint(Color.black);
//		plot.setOutlinePaint(Color.black);
//		plot.setDomainGridlinePaint(Color.white);
//		plot.setRangeGridlinePaint(Color.white);
	
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
		
//		XYSeriesCollection dataset1 = new XYSeriesCollection();
//		dataset1.addSeries(point);
//		plot.setDataset(0, dataset1);
//		XYSplineRenderer r1 = new XYSplineRenderer();
//		r1.setBaseLinesVisible(false);
//		plot.setRenderer(0, r1);
		
		XYSeriesCollection dataset2 = new XYSeriesCollection();
		XYSeries diag = new XYSeries(NP);
		diag.add(-2, -2);
		diag.add(26, 26);
		dataset2.addSeries(diag);
		
		XYSeriesCollection dataset3 = new XYSeriesCollection();
		XYSeries erroPos = new XYSeries("20%");
		erroPos.add(-2, -2*1.2);
		erroPos.add(26, 26*1.2);
		dataset3.addSeries(erroPos);
		
		XYSeriesCollection dataset4 = new XYSeriesCollection();
		XYSeries erroNeg = new XYSeries("20%");
		erroNeg.add(-2, -2*0.8);
		erroNeg.add(26, 26*0.8);
		dataset4.addSeries(erroNeg);
		
		plot.setDataset(i+1, dataset2);
		plot.setDataset(i+2, dataset3);
		plot.setDataset(i+3, dataset4);
		XYSplineRenderer r2 = new XYSplineRenderer();
		r2.setBaseShapesVisible(false);
		plot.setRenderer(i+1, r2);
		XYLineAndShapeRenderer r3 = new XYLineAndShapeRenderer();
		r3.setBaseShapesVisible(false);
		r3.setStroke(new BasicStroke(1.0f, BasicStroke.CAP_ROUND, BasicStroke.JOIN_ROUND,  
					1.0f, new float[] { 10.0f, 6.0f }, 0.0f));
		plot.setRenderer(i+2, r3);
		plot.setRenderer(i+3, r3);
				
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
}
