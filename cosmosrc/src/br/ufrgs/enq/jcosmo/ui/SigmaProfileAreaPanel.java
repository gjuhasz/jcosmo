/*
 * Copyright 2008 Rafael de Pelegrini Soares and Renan Pereira Gerber
 * 
 * This file is part of JCosmo.
 * 
 * JCosmo is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * JCosmo is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 * 
 * You should have received a copy of the GNU Lesser General Public License
 * along with JCosmo.  If not, see <http://www.gnu.org/licenses/>.
 */

package br.ufrgs.enq.jcosmo.ui;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Font;

import javax.swing.JPanel;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.data.Range;
import org.jfree.data.xy.DefaultTableXYDataset;
import org.jfree.data.xy.XYSeries;
import org.jfree.ui.RectangleInsets;

/**
 * Panel for sigma profile charts.
 * 
 * <p>This will create a stacked area plot for pure substances.
 * The {@link #addProfile(String, double[], double[])} adds a new descriptor for
 * the sigma profile. The multiple descriptors added are then stacked to form the
 * sigma profile.
 *  
 * @author Rafael de Pelegrini Soares
 */
public class SigmaProfileAreaPanel extends JPanel {
	private static final long serialVersionUID = 1L;

	XYPlot sigmaProfilePlot;
	JFreeChart sigmaProfileChart;
	DefaultTableXYDataset dataset = new DefaultTableXYDataset();
	
	public SigmaProfileAreaPanel() {
		setLayout(new BorderLayout());

		sigmaProfileChart = ChartFactory.createStackedXYAreaChart(null, 
				"sigma [e/Ų]", "p(sigma)�area [Ų]", dataset, PlotOrientation.VERTICAL, true, true, false);
		
		sigmaProfilePlot = sigmaProfileChart.getXYPlot();
		sigmaProfilePlot.getDomainAxis().setAutoRange(false);
		sigmaProfilePlot.getDomainAxis().setRange(new Range(-0.025, 0.025));
		sigmaProfilePlot.setBackgroundPaint(Color.WHITE);
		
		sigmaProfilePlot.setAxisOffset(new RectangleInsets(8, 0, 0, 16));
		sigmaProfilePlot.setDomainGridlinesVisible(false);
		sigmaProfilePlot.setRangeGridlinesVisible(false);
		
		Font font = new Font("DeJaVu Serif", Font.BOLD, 16);
		Font fontMini = new Font("DeJaVu Serif", 0, 12);
		
		sigmaProfilePlot.getDomainAxis().setLabelFont(font);
		sigmaProfilePlot.getRangeAxis().setLabelFont(font);
		sigmaProfilePlot.getDomainAxis().setTickLabelFont(fontMini);
		sigmaProfilePlot.getRangeAxis().setTickLabelFont(fontMini);
		
		sigmaProfilePlot.getDomainAxis().setAxisLinePaint(Color.BLACK);
		sigmaProfilePlot.getDomainAxis().setTickMarkPaint(Color.BLACK);
		sigmaProfilePlot.getDomainAxis().setTickMarkInsideLength(4);
		sigmaProfilePlot.getDomainAxis().setTickMarkOutsideLength(0);

		sigmaProfilePlot.getRangeAxis().setAxisLinePaint(Color.BLACK);
		sigmaProfilePlot.getRangeAxis().setTickMarkPaint(Color.BLACK);
		sigmaProfilePlot.getRangeAxis().setTickMarkInsideLength(4);
		sigmaProfilePlot.getRangeAxis().setTickMarkOutsideLength(0);
		
		
		add(new ChartPanel(sigmaProfileChart), BorderLayout.CENTER);
	}
	
	/**
	 * Remove all current profiles previouslly added.
	 */
	public void clearProfiles(){
		dataset.removeAllSeries();
	}
	
	/**
	 * Adds a new descriptor profile
	 * @param label the descriptor label
	 * @param sigma the sigma (ordered)
	 * @param area the area for each sigma
	 */
	public void addProfile(String label, double[] sigma, double[] area){
		int n = sigma.length;
		XYSeries comp= new XYSeries(label, false, false);

		// charges represent the center of the segments
		comp.add(sigma[0], area[0]);
		for(int j=1; j<n; ++j){
			double dsigma = (sigma[j]-sigma[j-1])/2;
			comp.add(sigma[j]-dsigma, area[j]);
			comp.add(sigma[j]+dsigma+1e-6, area[j]);
		}
		dataset.addSeries(comp);
		int series = dataset.getSeriesCount()-1;
		int start = 80, delta = 40;
		int rgb = start + delta*series;
		sigmaProfilePlot.getRenderer().setSeriesPaint(series, new Color(rgb, rgb, rgb));
	}
}
