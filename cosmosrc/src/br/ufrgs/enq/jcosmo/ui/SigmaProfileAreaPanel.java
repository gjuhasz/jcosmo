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
 * Dialog for building charts of the activity coefficient using COSMO-SAC model.
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
		
		sigmaProfilePlot.setAxisOffset(new RectangleInsets(0, 0, 0, 16));
		sigmaProfilePlot.getDomainAxis().setAxisLinePaint(Color.BLACK);
		sigmaProfilePlot.getRangeAxis().setAxisLinePaint(Color.BLACK);
		sigmaProfilePlot.getDomainAxis().setTickMarkPaint(Color.BLACK);
		sigmaProfilePlot.getRangeAxis().setTickMarkPaint(Color.BLACK);
		
		add(new ChartPanel(sigmaProfileChart), BorderLayout.CENTER);
	}
	
	public void clearProfiles(){
		dataset.removeAllSeries();
	}
	
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
		int start = 80, delta = 50;
		int rgb = start + delta*series;
		sigmaProfilePlot.getRenderer().setSeriesPaint(series, new Color(rgb, rgb, rgb));
	}
}