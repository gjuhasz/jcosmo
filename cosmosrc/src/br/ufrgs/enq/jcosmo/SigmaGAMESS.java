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

package br.ufrgs.enq.jcosmo;

import java.util.Scanner;

/**
 * Application for printing the sigma profile from a GAMESS LOG file.
 * 
 * @author rafael
 *
 */
public class SigmaGAMESS {
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		String filename;

		try {
			int nseg = 51;
			if(args.length > 0){
				filename = args[0];
				if(args.length>1){
					try{
						nseg = Integer.parseInt(args[1]);
					}
					catch (NumberFormatException e) {
						System.out.println("Invalid number of segments");
						return;
					}
				}			
			}
			else{
				Scanner in = new Scanner(System.in);
				System.out.print("Enter with the GAMESS LOG FILE: ");
				filename = in.nextLine();
			}
			SigmaProfileGenerator a = new SigmaProfileGenerator(SigmaProfileGenerator.FileType.GAMESS, nseg);
			a.parseFile(filename);
			a.printProfile(System.out);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
