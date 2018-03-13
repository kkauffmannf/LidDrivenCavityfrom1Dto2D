/*
 * readinput.cpp
 *
 * Reads a file named input.txt
 * Its format is
 *
 * var1 = value
 * var2 = value
 * .
 * .
 * .
 *
 *
 *  Created on: May 23, 2017
 *      Author: Karla Kauffmann
 *      E-mail: karla.kauffmann@cchen.cl
 */

#include <iostream>
#include <string>
#include <fstream>
#include <sstream>

#include "globals.h"

using namespace std;

void readinput()
{
	   /* Opens, reads and stores values from input.txt file */
		   ifstream inFile("input.txt");
		   string strOneLine;
		   string temp;
		   string delimiter = " = "; /* reads value after the equal sign */
		   size_t pos = 0;
		   string token;
		   int counter=0;

		   while (inFile)
		   {
		      getline(inFile, strOneLine); /* reads a line */

		      while ((pos = strOneLine.find(delimiter)) != string::npos) {

		    	 /* deletes the text in line that goes before the equal sign */
		         strOneLine.erase(0, pos + delimiter.length());

		         /* stores the value after the equal sign */
		         temp = strOneLine;

		      }

		      /* Stores the values from input.txt in the corrresponding variables */
		      if (counter == 0) {
		        istringstream ( temp ) >> Nodesx;
		      }
		      if (counter == 1) {
		        istringstream ( temp ) >> Nodesy;
		      }
		      if (counter == 2) {
		        istringstream ( temp ) >> Lx;
		      }
		      if (counter == 3) {
		        istringstream ( temp ) >> Ly;
		      }
		      if (counter == 4) {
		        istringstream ( temp ) >> p_init;
		      }
		      if (counter == 5) {
		        istringstream ( temp ) >> density;
		      }
		      if (counter == 6) {
		      	istringstream ( temp ) >> lid_velocity;
		      }
		      if (counter == 7) {
		    	istringstream ( temp ) >> Reynolds_num;
		      }
		      if (counter == 8) {
		    	istringstream ( temp ) >> urfu;
		      }
		      if (counter == 9) {
		      	istringstream ( temp ) >> urfv;
		      }
		      if (counter == 10) {
		    	istringstream ( temp ) >> urfp;
		      }
		      if (counter == 11) {
		    	istringstream ( temp ) >> residual_threshold;
		      }
		      counter++;
		   }


		   inFile.close();

	       /* prints out the input values */
		   cout << "These are the input values:" << endl;
		   cout << "Nx = " << Nodesx << endl;
		   cout << "Ny = " << Nodesy << endl;
		   cout << "Lx = " << Lx << endl;
		   cout << "Ly = " << Ly << endl;
		   cout << "p_init = " << p_init << endl;
		   cout << "density = " << density << endl;
		   cout << "lid_velocity = " << lid_velocity << endl;
		   cout << "Reynolds_num = " << Reynolds_num << endl;
	       cout << "urfu = " << urfu << endl;
	       cout << "urfv = " << urfv << endl;
		   cout << "urfp = " << urfp << endl;
		   cout << "residual_threshold = "  << residual_threshold << endl;
		   cout << endl;

}
