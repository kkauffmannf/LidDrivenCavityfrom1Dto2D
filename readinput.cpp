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
	        istringstream ( temp ) >> N;
	      }
	      if (counter == 1) {
	        istringstream ( temp ) >> L;
	      }
	      if (counter == 2) {
	        istringstream ( temp ) >> p_in;
	      }
	      if (counter == 3) {
	        istringstream ( temp ) >> p_out;
	      }
	      if (counter == 4) {
	        istringstream ( temp ) >> density;
	      }
	      if (counter == 5) {
	        istringstream ( temp ) >> Area_in;
	      }
	      if (counter == 6) {
	      	istringstream ( temp ) >> Area_out;
	      }
	      if (counter == 7) {
	    	istringstream ( temp ) >> m_dot;
	      }
	      if (counter == 8) {
	    	istringstream ( temp ) >> urfu;
	      }
	      if (counter == 9) {
	    	istringstream ( temp ) >> urfp;
	      }
	      if (counter == 10) {
	    	istringstream ( temp ) >> residual_threshold;
	      }
	      counter++;
	   }


	   inFile.close();

       /* prints out the input values */
	   cout << "These are the input values:" << endl;
	   cout << "N = " << N << endl;
	   cout << "L = " << L << endl;
	   cout << "p_in = " << p_in << endl;
	   cout << "p_out = " << p_out << endl;
	   cout << "density = " << density << endl;
	   cout << "Area_in = " << Area_in << endl;
	   cout << "Area_out = " << Area_out << endl;
       cout << "m_dot = " << m_dot << endl;
	   cout << "urfu = " << urfu << endl;
	   cout << "urfp = " << urfp << endl;
	   cout << "residual_threshold = " << residual_threshold << endl;
	   cout << endl;

}
