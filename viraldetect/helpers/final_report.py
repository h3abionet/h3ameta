#!/usr/bin/env python2.7
import os
import sys 
import re
import glob
import csv
from pandas import DataFrame
import pandas as pd
import numpy as np
import webbrowser
import pandas



def main(argv):
    if len(argv) >= 2:
        print ("\nUsage: python2.7 final_report.py <krakenfile> <mappingstatoutput1> <mappingstatoutput2> <mappingstatoutput3> <kronaoutput> <outputfile>")
        sys.exit(2)

kraken = sys.argv[1]
statH =  sys.argv[2]+ "bammappingStats.txt"
statE =  sys.argv[3]+ "bammappingStats.txt"
statC =  sys.argv[4]+ "bammappingStats.txt"
krona =  sys.argv[5]
outputfileconcat = sys.argv[6]



stringToMatch = ['hepatitis', 'Ebola' , 'Chikungunya' ]
df = []
for key in stringToMatch :
	#print (key)
	krakenfile = open (kraken , "r")	
	trouv = key,"Not Found"
	for line in krakenfile :	
		if key in line: 
			trouv = key,"Found"
	#print (trouv)

	df.append(trouv)
#print(df)
ff = DataFrame.from_records(df) 
#print (ff)
#bwa output
import glob

bwa = glob.glob("*bammappingStats.txt")
#print (bwa)
search = ['mapped (']
df2 = []
for file in bwa :
    with open(file, 'r') as f:
        ligne = f.readlines()
        #print (ligne)
        for line in ligne:
            if "mapped (" in line :
                match = (re.sub('\.bammappingStats.txt$', '', file),line)
                #print (match)
                df2.append(match)
                #print(df2)
gg = DataFrame.from_records(df2)
result = pd.concat([gg, ff], axis=1, ignore_index=True)
result = result.drop([2], axis = 1)
result.columns = ['Pathogen', '%mapped' , 'kraken']
#print(result)

result.to_html("Table.html")

#convert kraken tsv to html :#https://github.com/dbohdan/csv2html?fbclid=IwAR0L3a41xXXVpWXxH-q-httRuB1x-6S6KXQysBAx5dGhBaE3Ntg6WkVs7yc
table = ""
nargs = len(sys.argv)
with open(kraken, 'rb') as tsvfile :
	csv_table=pd.read_table(tsvfile , sep='\t')
	csv_table.columns =['Percentage of reads covered by the clade rooted at this taxon' ,'Number of reads covered by the clade rooted at this taxon','Number of reads assigned directly to this taxon','Rank code','NCBI taxonomy ID','scientific name']	
 #convert tsv to csv first
	
	csv_table.to_csv('kraken.csv.tmp' )	
	#print  (csv_table)
	table = ""
	with open('kraken.csv.tmp', 'rb') as csvfile:
		
		reader = csv.reader(csvfile , delimiter=',')
		#print (reader)

		table += "<table border='1'>\n"
		for row in reader:
			table += "<tr>\n" + "".join(["<td>%s</td>\n" % 
				     item for item in row]) + "</tr>\n"
		table += "</table>\n"

	# convert csv to html
	if nargs > 2:
		with open(sys.argv[6], 'w') as htmlfile:
			htmlfile.write(table)
			#print type (htmlfile)
		
	else:
		output_file.write(table)
#                                                                                                                                                                         print (table)

#URL =   "input file name is " % kraken
#print (URL)
#message = ""
#new_message = message.format(URL="Kraken")

######################"report
#f = open ('final_report.html', "w")


with open(outputfileconcat , "w") as f :
	message1 = """<html>
    
    <head>
        <header>
            <h1 align="center"> General report</h1>
	 <h1 align="center"> <p> Outputfile name : """
	
	message2 = outputfileconcat
	
	message3 = """ </h1></p>
        </header>
    </head>
    <head>
        <title> Input file name </title>
	
    </head>
    <body><p>
        <section id="krona">
            <h2 id> Interactive vizualisation of taxa relative abundance by Krona </h2>
                <p>Krona result are displayed here :</p> <a href ='krona.html'> krona </a>
		
        </section></p></body>
        </section>
    
    </section>
    <h2> Best Hits </h2>
        <section id="hits">
            <embed src='Table.html' width="800px" height="300px" />
            <div class="sec">
            <div class="leftexp">
            </div>
    </section>
    
    <section id="kraken">
        <h2> kraken_report </h2>
            <embed src='kraken_report.html' width="1600px" height="1000px" />
            <div class="sec">
            <div class="leftexp">
            </div>
    </section>
    </body>
    </html>"""
	f.write(message1 + message2+ message3)
	f.close()

filenname =  outputfileconcat

webbrowser.open_new_tab(filenname)


