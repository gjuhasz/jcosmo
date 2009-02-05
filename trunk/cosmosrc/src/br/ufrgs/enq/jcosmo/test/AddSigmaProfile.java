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
public class AddSigmaProfile {

	private static final String DRIVER = "org.apache.derby.jdbc.EmbeddedDriver";

	public static void main(String[] args) {
		try {

			if(args.length!=5){
				System.out.println("input argumets are: NAME FORMULA CAS SIGMA_FILE CAV_VOLUME");
				return;
			}
			String name = args[0];
			String formula = args[1];
			String CAS = args[2];
			String sigmaFile = args[3];
			String cavVolume = args[4];

			String dbName = "../cosmodb";
			//		String password = "br.eng.rps";

			String connectionURL = "jdbc:derby:" + dbName;
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

			addProfile(statement, formula, name, CAS, sigmaFile, cavVolume);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	static void addProfile(Statement statement, String formula, String name, String CAS, String sigmaFile,
			String cavVolume) throws Exception{

		ArrayList<Double> charge = new ArrayList<Double>();
		ArrayList<Double> sigma = new ArrayList<Double>();
		readSigmaProfile(sigmaFile, charge, sigma);

		int id = -1;
		String count = "SELECT id from cosmoindex where name like " + quote(name);
		ResultSet res = statement.executeQuery(count);
		if(res.next()){
			System.out.println(name + " already in database");
			return;
		}

		count = "SELECT count(id) from cosmoindex";
		res = statement.executeQuery(count);
		res.next();
		id = res.getInt(1)+1;

		double dcavVolume = Double.parseDouble(cavVolume)*0.529177249*0.529177249*0.529177249;

		String sql = "INSERT INTO COSMOINDEX (ID, FORMULA, NAME, CAS, Vcosmo) VALUES ("+
		id + ',' + quote(formula) + ',' + quote(name) + ',' + quote(CAS) + ',' + dcavVolume + ")";

		statement.execute(sql);

		count = "SELECT count(id) from SIGMAPROFILE";
		res = statement.executeQuery(count);
		res.next();
		int sid = res.getInt(1)+1;
		for (int j = 0; j < charge.size(); j++) {
			sql = "INSERT INTO SIGMAPROFILE (ID, CompID, CHARGE, SIGMA) VALUES(" +
			sid++ + ',' + id + ',' + charge.get(j) + ',' + sigma.get(j) + ")";
			statement.execute(sql);
		}
	}

	public static String quote(String in){
		if(in.length() == 0)
			return "NULL";
		if(in.startsWith("\""))
			return in.replace("\"", "'");
		else
			return "'" + in + "'";
	}

	public static void readSigmaProfile(String sigmaFile, ArrayList<Double> charge, ArrayList<Double> sigma) throws Exception{
		Scanner input = new Scanner(new File(sigmaFile));
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
