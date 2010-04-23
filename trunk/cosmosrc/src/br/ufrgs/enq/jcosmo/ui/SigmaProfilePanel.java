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

import java.awt.BasicStroke;
import java.awt.BorderLayout;

import javax.swing.JPanel;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYStepRenderer;
import org.jfree.data.Range;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

/**
 * Dialog for building charts of the activity coefficient using COSMO-SAC model.
 *  
 * @author Rafael de Pelegrini Soares
 */
public class SigmaProfilePanel extends JPanel {
	private static final long serialVersionUID = 1L;

	XYPlot sigmaProfilePlot;
	JFreeChart sigmaProfileChart;
	XYSeriesCollection dataset = new XYSeriesCollection();
	
	XYStepRenderer stepRenderer = new XYStepRenderer();
	
	boolean err = false;

	private BasicStroke dashed;

	public SigmaProfilePanel(String title) {
		setLayout(new BorderLayout());

		sigmaProfileChart = ChartFactory.createXYLineChart(title, 
				"sigma", "P^x", dataset, PlotOrientation.VERTICAL, true, true, false);
		sigmaProfilePlot = sigmaProfileChart.getXYPlot();
		sigmaProfilePlot.getDomainAxis().setAutoRange(false);
		sigmaProfilePlot.getDomainAxis().setRange(new Range(-0.025, 0.025));
		
		add(new ChartPanel(sigmaProfileChart), BorderLayout.CENTER);
		
		dashed = new BasicStroke(1.0f, BasicStroke.CAP_ROUND, BasicStroke.JOIN_ROUND,  
				1.0f, new float[] { 10.0f, 6.0f }, 0.0f);
		
//		sigmaProfilePlot.setBackgroundPaint(Color.lightGray);
//		sigmaProfilePlot.setDomainGridlinePaint(Color.white);
//		sigmaProfilePlot.setRangeGridlinePaint(Color.white);
	}
	
	public void clearProfiles(){
		dataset.removeAllSeries();
	}
	
	public void addProfile(String label, double[] sigma, double[] area){
		addProfile(label, sigma, area, false);
	}
	

	public void addProfile(String label, double[] sigma, double[] area, boolean dashed){
		int n = sigma.length;
		XYSeries comp= new XYSeries(label);

		// charges represent the center of the segments
		comp.add(sigma[0], area[0]);
		for(int j=1; j<n; ++j){
			comp.add(sigma[j]-(sigma[j]-sigma[j-1])/2, area[j]);
		}
		dataset.addSeries(comp);
		sigmaProfilePlot.setRenderer(dataset.getSeriesCount()-1, stepRenderer);
		if(dashed)
			sigmaProfilePlot.getRenderer().setSeriesStroke(dataset.getSeriesCount()-1, this.dashed);
		else
			sigmaProfilePlot.getRenderer().setSeriesStroke(dataset.getSeriesCount()-1, new BasicStroke(2.5f));
	}
}
