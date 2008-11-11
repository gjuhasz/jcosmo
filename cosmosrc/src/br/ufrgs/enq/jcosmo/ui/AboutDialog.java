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
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.JButton;
import javax.swing.JDialog;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;

/**
 * The About dialog.
 * 
 * @author rafael
 *
 */
public class AboutDialog extends JDialog {
	private static final long serialVersionUID = 1L;

	public AboutDialog(JFrame parent) {
		super(parent, "About this program...", true);

		JLabel text = new JLabel(
				"<html><center>This is a small program for the calculation of" +
				"<br>activity coefficients and g<sup>E</sup> of binary mixtures" +
				"<br>using a COSMO-SAC implementation which is" +
				"<br>a re-implementation of COSMO-RS by A. Klamt." +
				
				"<br><br>For more details check http://code.google.com/p/jcosmo" +
				
				"<br><br>This software is provided AS IS" +
				"<br>WITHOUT WARRANTY OF ANY KIND (LGPL license)." + 

				"<br><br>This implementation is inspired on the code" +
				"<br>available at http://www.design.che.vt.edu." +
				
				"<br><br>Rafael de Pelegrini Soares - http://www.rps.eng.br" +
				"<br>Renan Pereira Gerber" +
				"<br> $Date: 2008-03-29 23:49:38 -0300 (Sat, 29 Mar 2008) $" +
				
				"</center></html>");
		getContentPane().add(text, BorderLayout.CENTER);

		JPanel p2 = new JPanel();
		JButton ok = new JButton("Ok");
		p2.add(ok);
		getContentPane().add(p2, BorderLayout.SOUTH);

		ok.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent evt) {
				setVisible(false);
			}
		});

		pack();
		setLocationRelativeTo(parent);
		setVisible(true);
	}
}
