package br.ufrgs.enq.jcosmo2;

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

import br.ufrgs.enq.jcosmo2.SigmaProfileGenerator2.FileType;

/**
 * Diagonal for averaging energy effects.
 * 
 * @author rafael
 *
 */
public class AveragingEnergyCorrected3 extends JFrame implements XYToolTipGenerator{
	private static final long serialVersionUID = 1L;

	List<String> tps = new Vector<String>();

	List<Color> color = new ArrayList<Color>();

	public AveragingEnergyCorrected3() {

		// Configuration
		String folder = "profiles/RM1_1.18/";
		String extension = ".cos";
		FileType type = SigmaProfileGenerator2.FileType.MOPAC;
		double rav = 1.5;

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
		
		XYSeries points0 = new XYSeries("Points", false);
		dataset.addSeries(points0);
		XYSeries points1 = new XYSeries("Points", false);
		dataset.addSeries(points1);
		r.setSeriesLinesVisible(0, false);
		r.setSeriesShapesVisible(0, true);
		r.setSeriesToolTipGenerator(0, this);
		r.setSeriesLinesVisible(1, false);
		r.setSeriesShapesVisible(1, true);
		r.setSeriesToolTipGenerator(1, this);
		
		ChartPanel chartPanel = new ChartPanel(chart);
		JPanel jp = new JPanel(new BorderLayout());
		jp.add(chartPanel);
		
		add(jp, BorderLayout.CENTER);
		setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		
		System.out.println("Name\tEnergy\tEnergy_avg");
		
		SimpleRegression line = new SimpleRegression();
		SimpleRegression line0 = new SimpleRegression();
		
		File folderObj = new File(folder);
		String names[] = folderObj.list();
		
		double maxEnergy = -Double.MAX_VALUE;
		for(String name : names){
			if(!name.endsWith(extension))
				continue;

			name = name.substring(0, name.length() - extension.length());
			String fileName = folder + name + extension;

			SigmaProfileGenerator2 sparser = new SigmaProfileGenerator2(type);
			try {
				sparser.parseFile(fileName, rav);
			} catch (Exception e) {
				e.printStackTrace();
				continue;
			}

			double []sigmaAvg = sparser.getAveragedChargeDensity();
			double []area = sparser.getOriginalArea();
			double []sigma = sparser.getOriginalChargeDensity();

			double energy = 0, energyAvgCorr = 0, energyAvg0 = 0;
			double areaT = 0;
			for (int i = 0; i < sigmaAvg.length; i++) {
				areaT += area[i];
				
				double sigmai = sigma[i];
				double sigmaAvgi = sigmaAvg[i];

				energy += area[i] * sigmai*sigmai;
				energyAvg0 += area[i] * sigmaAvgi*sigmaAvgi;
			}
			
			for (int i = 0; i < sigmaAvg.length; i++) {
				//put back the energy as a factor
				double sigmaAvgCorri = sigmaAvg[i]*Math.sqrt(energy/energyAvg0);
				
				energyAvgCorr += area[i] * sigmaAvgCorri*sigmaAvgCorri;
			}

			points0.add(energy/areaT*1e3, energyAvgCorr/areaT*1e3);
			points1.add(energy/areaT*1e3, energyAvg0/areaT*1e3);
			line.addData(energy/areaT*1e3, energyAvgCorr/areaT*1e3);
			line0.addData(energy/areaT*1e3, energyAvg0/areaT*1e3);
			tps.add(name + ", " + energy/areaT*1e3);
			
			maxEnergy = Math.max(maxEnergy, energy/areaT*1e3);

			System.out.println(name + "\t" + energy/areaT*1e3 + "\t" + energyAvgCorr/areaT*1e3);
		}

		
		XYSeries diag = new XYSeries("Slope = " + line.getSlope() + ", R2 = " + line.getRSquare() + ", MSE = " + line.getMeanSquareError());
		diag.add(0, line.predict(0));
		diag.add(maxEnergy*1.2, line.predict(maxEnergy*1.2));
		dataset.addSeries(diag);
		r.setSeriesStroke(2, dashed);
		r.setSeriesLinesVisible(2, true);
		
		XYSeries diag0 = new XYSeries("Slope = " + line0.getSlope() + ", R2 = " + line0.getRSquare() + ", MSE = " + line0.getMeanSquareError());
		diag0.add(0, line0.predict(0));
		diag0.add(maxEnergy*1.2, line0.predict(maxEnergy*1.2));
		dataset.addSeries(diag0);
		r.setSeriesStroke(3, dashed);
		r.setSeriesLinesVisible(3, true);
		
		
		setSize(800,600);
		setVisible(true);
	}

	public String generateToolTip(XYDataset dataset, int series, int item) {
		return tps.get(item);
	}
	
	public static void main(String[] args) {
		new AveragingEnergyCorrected3();
	}
}
