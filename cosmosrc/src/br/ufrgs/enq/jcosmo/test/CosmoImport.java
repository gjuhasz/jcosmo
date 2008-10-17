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

package br.ufrgs.enq.jcosmo.test;

import java.io.File;
import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.ResultSet;
import java.sql.Statement;
import java.util.ArrayList;
import java.util.Locale;
import java.util.Scanner;

import com.csvreader.CsvReader;

/**
 * Class used to import data from COSMO VT2005 database (pure text files)
 * into a SQL database.
 * 
 * <p>More preciselly an Apache Derby database is created. This code should
 * not be executed again if the import was already done (the code is kept
 * only for historical purposes).
 * 
 * @author rafael
 *
 */
public class CosmoImport {

	private static final String DRIVER = "org.apache.derby.jdbc.EmbeddedDriver";

	public static void main(String[] args) {
		CosmoImport cosmoImport = new CosmoImport();

		try {
			cosmoImport.importDababase();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	private String quote(String in){
		if(in.length() == 0)
			return "NULL";
		if(in.startsWith("\""))
			return in.replace("\"", "'");
		else
			return "'" + in + "'";
	}

	public void importDababase() throws Exception{

		String dbName = "cosmodb";
//		String password = "br.eng.rps";

		String connectionURL = "jdbc:derby:" + dbName + ";create=true";
//		+ ";dataEncryption=true;bootPassword=" + password;

		/*
		 **  Load the Derby driver. 
		 **     When the embedded Driver is used this action start the Derby engine.
		 **  Catch an error and suggest a CLASSPATH problem
		 */
		Class.forName(DRIVER); 
		System.out.println(DRIVER + " loaded. ");

		// Create (if needed) and connect to the database
		Connection conn = DriverManager.getConnection(connectionURL);		 
		System.out.println("Connected to database " + dbName);

		Statement statement = conn.createStatement(ResultSet.TYPE_SCROLL_INSENSITIVE, ResultSet.CONCUR_READ_ONLY);

		CsvReader reader = new CsvReader("COSMO/VT-2005_Sigma_Profile_Database_Index_v2.csv");
		reader.readHeaders();

		int id = 0, segments;
		String name, formula, CAS, family, preOptTool;
		double Vcosmo, T, lnPvap;
		ArrayList<Double> charge = new ArrayList<Double>();
		ArrayList<Double> sigma = new ArrayList<Double>();

		String sql = "CREATE TABLE COSMOINDEX (ID INTEGER PRIMARY KEY,Formula VARCHAR(128),Name VARCHAR(128),CAS VARCHAR(128), Family VARCHAR(128)" +
		", Vcosmo DOUBLE,Segments INTEGER,preOptTool VARCHAR(128),T DOUBLE,lnPvap DOUBLE)";
		statement.execute(sql);

		sql = "CREATE TABLE SIGMAPROFILE (ID INTEGER PRIMARY KEY, CompID INTEGER, Charge DOUBLE, Sigma DOUBLE)";
		statement.execute(sql);

		int sigmaid = 0;

		while(reader.readRecord()){
			String[] values = reader.getValues();

			int i = 0;
			id = Integer.parseInt(values[i++]);
			formula = quote(values[i++]);
			name = quote(values[i++]);
			CAS = quote(values[i++]);
			family = quote(values[i++]);
			Vcosmo = Double.parseDouble(values[i++]);
			segments = Integer.parseInt(values[i++]);
			preOptTool = quote(values[i++]);
			T = Double.parseDouble(values[i++]);
			lnPvap = Double.parseDouble(values[i++]);
			sql = "INSERT INTO COSMOINDEX VALUES ("+
			id + ',' + formula + ',' + name + ',' + CAS +
			',' + 	family + ',' + Vcosmo + ',' + segments + ',' + preOptTool +
			',' + T + ',' + lnPvap +
			")";
			
			System.out.println("Adding compound:" + name);

			statement.execute(sql);

			readSigmaProfile(id, charge, sigma);

			for (int j = 0; j < charge.size(); j++) {
				sql = "INSERT INTO SIGMAPROFILE VALUES(" +
				sigmaid++ + ',' + id + ',' +
				charge.get(j) + ',' + sigma.get(j) + ")";
					statement.execute(sql);
			}
		}
		reader.close();
		
//		sql = "SELECT COUNT(ID) FROM COSMOINDEX";
//		statement.execute(sql);
//		ResultSet result = statement.getResultSet();
//		System.out.println("Components added:" + result.getInt(0));
//		
//		sql = "SELECT COUNT(ID) FROM SIGMAPROFILE";
//		statement.execute(sql);
//		result = statement.getResultSet();
//		System.out.println("Sigma points added:" + result.getInt(0));
	}

	public static void readSigmaProfile(int id, ArrayList<Double> charge, ArrayList<Double> sigma) throws Exception{
		String filename = String.format("COSMO/VT2005/VT2005-%4d-PROF.txt", id).replace(' ', '0');

		Scanner input = new Scanner(new File(filename));
		input.useLocale(Locale.US);

		charge.clear();
		sigma.clear();

		double chargei, sigmai;
		while(input.hasNext()){
			chargei = input.nextDouble();
			sigmai = input.nextDouble();

			charge.add(chargei);
			sigma.add(sigmai);
		}
		input.close();
	}
}
