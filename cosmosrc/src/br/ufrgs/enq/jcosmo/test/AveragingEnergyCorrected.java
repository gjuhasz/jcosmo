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

import br.ufrgs.enq.jcosmo.SigmaProfileGenerator;
import br.ufrgs.enq.jcosmo.SigmaProfileGenerator.FileType;

/**
 * Diagonal for averaging energy effects.
 * 
 * @author rafael
 *
 */
public class AveragingEnergyCorrected extends JFrame implements XYToolTipGenerator{
	private static final long serialVersionUID = 1L;

	List<String> tps = new Vector<String>();

	List<Color> color = new ArrayList<Color>();

	public AveragingEnergyCorrected() {

		// Configuration
//		String folder = "profiles/RM1_1.18/";
//		String extension = ".cos";
//		FileType type = SigmaProfileGenerator.FileType.MOPAC;
//		double rav = 1.0;
//		double f_ortho = 0.8033901717857764;
//		double rav2 = 1.5*rav;
//		double f_corr = -2.0;
		
//		String folder = "profiles/RM1_1.18/";
//		String extension = ".cos";
//		FileType type = SigmaProfileGenerator.FileType.MOPAC;
//		double rav = 1.5;
//		double f_ortho = 0.6665617963244067;
//		double rav2 = 1.5*rav;
//		double f_corr = -2.0;
		
		String folder = "profiles/RM1_1.18/";
		String extension = ".cos";
		FileType type = SigmaProfileGenerator.FileType.MOPAC;
		double rav = 1.5;
		double f_ortho = 0.41830198448203804;
		double rav2 = 2*rav;
		double f_corr = -2.6;
		
//		String folder = "profiles/gamess/";
//		String extension = ".gout";
//		FileType type = SigmaProfileGenerator.FileType.GAMESS;
//		double rav = 1.5;
//		double f_ortho = 0.7241079673777269;
//		double rav2 = 1.5*rav;
//		double f_corr = -4.0;

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

			name = name.substring(0, name.length() - extension.length());
			String fileName = folder + name + extension;

			SigmaProfileGenerator sigmaParser;

			sigmaParser = new SigmaProfileGenerator(type);
			double[] sigmaAvg;
			double[] sigmaAvg2;
			try {
				sigmaParser.parseFile(fileName, rav);
				sigmaAvg = sigmaParser.getAveragedChargeDensity();
				
				sigmaParser.parseFile(fileName, rav2);
				sigmaAvg2 = sigmaParser.getAveragedChargeDensity();
			} catch (Exception e) {
				continue;
			}

			double[] area = sigmaParser.getOriginalArea();
			double []sigma = sigmaParser.getOriginalChargeDensity();


			double energy = 0, energyAvg = 0;
			double areaT = 0;
			for (int i = 0; i < sigmaAvg.length; i++) {
				areaT += area[i];
				
				// Orthogonalize the sigmaAvg2
				sigmaAvg2[i] = sigmaAvg2[i] - f_ortho*sigmaAvg[i];
				
				energy += area[i]*sigma[i]*sigma[i];
				energyAvg += area[i]*sigmaAvg[i]*(sigmaAvg[i] + f_corr*sigmaAvg2[i]);
			}
			points.add(energy/areaT*1e3, energyAvg/areaT*1e3);
			line.addData(energy/areaT*1e3, energyAvg/areaT*1e3);
			tps.add(name + ", " + energy*1e3);
			
			maxEnergy = Math.max(maxEnergy, energy/areaT*1e3);

			System.out.println(name + "\t" + energy/areaT + "\t" + energyAvg/areaT);
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
		new AveragingEnergyCorrected();
	}
}
