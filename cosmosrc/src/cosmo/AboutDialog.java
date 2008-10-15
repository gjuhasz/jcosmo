package cosmo;

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
				"<br>using the COSMO-SAC model." +
				
				"<br><br>This software is provided AS IS and should be" +
				"<br>used for evaluation purposes only." + 
				
				"<br><br>Rafael de Pelegrini Soares - http://www.rps.eng.br" +
				"<br>Renan Pereira Gerber" +
				"<br> $Date: 2008-03-29 23:49:38 -0300 (Sat, 29 Mar 2008) $" + 
				
				"<br><br>Based on code and data available at" +
				"<br>Virginia Tech - http://www.design.che.vt.edu" + 
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
