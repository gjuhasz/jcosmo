package br.ufrgs.enq.jcosmo.test;

import java.awt.BasicStroke;
import java.awt.BorderLayout;
import java.awt.Color;
import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.Vector;

import javax.swing.JFrame;
import javax.swing.JPanel;

import org.apache.commons.math.stat.regression.SimpleRegression;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.labels.XYToolTipGenerator;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import org.jfree.data.xy.XYDataset;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

import br.ufrgs.enq.jcosmo.COSMOSAC;
import br.ufrgs.enq.jcosmo.SigmaProfileGenerator;

/**
 * Diagonal for averaging energy effects.
 * 
 * @author rafael
 *
 */
public class AveragingEnergy extends JFrame implements XYToolTipGenerator{
	private static final long serialVersionUID = 1L;

	List<String> tps = new Vector<String>();

	List<Color> color = new ArrayList<Color>();

	public AveragingEnergy() {

		// Configuration
		String folder = "profiles/RM1/";
		String extension = ".cos";

		XYSeriesCollection dataset = new XYSeriesCollection();

		XYPlot plot;
		JFreeChart chart = ChartFactory.createXYLineChart(null, 
				"Energy", "Energy_avg", dataset,
				PlotOrientation.VERTICAL, true, true, false);
		plot = (XYPlot) chart.getPlot();		
//		plot.setBackgroundPaint(Color.black);
	
		XYLineAndShapeRenderer r = (XYLineAndShapeRenderer) plot.getRenderer();
		
		BasicStroke dashed = new BasicStroke(1.0f, BasicStroke.CAP_ROUND, BasicStroke.JOIN_ROUND,  
				1.0f, new float[] { 10.0f, 6.0f }, 0.0f);
		
//		XYSeries erroPos = new XYSeries("-1 ln unit");
//		erroPos.add(-2, -1);
//		erroPos.add(26, 27);
//		dataset.addSeries(erroPos);
//		r.setSeriesStroke(i, dashed);
//		r.setSeriesLinesVisible(i++, true);
//		
//		XYSeries erroNeg = new XYSeries("+1 ln unit");
//		erroNeg.add(-2, -3);
//		erroNeg.add(26, 25);
//		dataset.addSeries(erroNeg);
//		r.setSeriesStroke(i, dashed);
//		r.setSeriesLinesVisible(i++, true);
		
		XYSeries points = new XYSeries("Points", false);
		dataset.addSeries(points);
		r.setSeriesLinesVisible(0, false);
		r.setSeriesShapesVisible(0, true);
		r.setSeriesToolTipGenerator(0, this);
		
		ChartPanel chartPanel = new ChartPanel(chart);
		JPanel jp = new JPanel(new BorderLayout());
		jp.add(chartPanel);
		
		add(jp, BorderLayout.CENTER);
		setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		
		System.out.println("Name\tEnergy\tEnergy_avg");
		
		SimpleRegression line = new SimpleRegression();
		
		File folderObj = new File(folder);
		String names[] = folderObj.list();
		
		double maxEnergy = -Double.MAX_VALUE;
		for(String name : names){
			if(!name.endsWith(extension))
				continue;

			name = name.substring(0, name.length() - 4);
			String fileName = folder + name + extension;

			SigmaProfileGenerator sigmaParser;

			sigmaParser = new SigmaProfileGenerator(SigmaProfileGenerator.FileType.MOPAC);
			try {
				sigmaParser.parseFile(fileName, COSMOSAC.RAV);
			} catch (Exception e) {
				continue;
			}

			double[] area = sigmaParser.getOriginalArea();
			double []sigma = sigmaParser.getOriginalChargeDensity();
			double[] sigmaAvg = sigmaParser.getAveragedChargeDensity();

			double energy = 0, energyAvg = 0;
			for (int i = 0; i < sigmaAvg.length; i++) {
				energy += area[i]*sigma[i]*sigma[i];
				energyAvg += area[i]*sigmaAvg[i]*sigmaAvg[i];
			}
			points.add(energy*1e3, energyAvg*1e3);
			line.addData(energy*1e3, energyAvg*1e3);
			tps.add(name + ", " + energy*1e3);
			
			maxEnergy = Math.max(maxEnergy, energy*1e3);

			System.out.println(name + "\t" + energy + "\t" + energyAvg);
		}

		
		XYSeries diag = new XYSeries("Slope = " + line.getSlope() + ", R2 = " + line.getR() + ", MSE = " + line.getMeanSquareError());
		diag.add(0, line.predict(0));
		diag.add(maxEnergy*1.2, line.predict(maxEnergy*1.2));
		dataset.addSeries(diag);
		r.setSeriesStroke(1, dashed);
		r.setSeriesLinesVisible(1, true);
		
		
		setSize(800,600);
		setVisible(true);	
	}

	public String generateToolTip(XYDataset dataset, int series, int item) {
		if(series == 0)
			return tps.get(item);
		return null;
	}
	
	public static void main(String[] args) {
		new AveragingEnergy();
	}
}
