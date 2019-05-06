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
import argparse


kraken = sys.argv[1]
stat =  sys.argv[2]
krona =  sys.argv[3]
outputfileconcat = sys.argv[4]

#==================================================================================================================
#Detecttion of viruses
stringToMatch = ['hepatitis', 'Ebola' , 'Chikungunya' ]
df = []
for key in stringToMatch :
    #print (key)
    krakenfile = open (kraken , "r")
    trouv = key,"Not Found"
    for line in krakenfile :
        if key in line:
            trouv = (key,"Found")
#print (trouv)
    df.append(trouv)
#print(df)
kr_df = DataFrame.from_records(df)
#print (kr_df)

#===================================================================================
#bwa output

searchHep = "FJ407092.1"
searchEbo= "KP096420.1"
searchChik = "MG049915.1"

stat_df = pd.read_csv(stat, sep="\t", header=None)
stat_df.columns = ["a", "b", "c", "d"]
stat_df= stat_df.drop([3], axis=0) #delete last line
#print (stat_df)

#merge results
result = pd.concat([kr_df, stat_df], axis=1)
#print(result)

result.columns = [ 'Pathogen', 'kraken result', 'Reference sequence identifier' ,'Reference sequence length', 'Number of mapped reads',  'Number of unmapped reads']


result.to_html("Table.html" , index=False)

#=========================================================================================================================
#convert kraken tsv to html :#https://github.com/dbohdan/csv2html?fbclid=IwAR0L3a41xXXVpWXxH-q-httRuB1x-6S6KXQysBAx5dGhBaE3Ntg6
#add header
table = ""
nargs = len(sys.argv)
with open(kraken, 'rb') as tsvfile :
        csv_table=pd.read_csv(tsvfile , sep='\t')
        csv_table.columns = ['Percentage of reads covered by the clade rooted at this taxon' ,'Number of reads covered by the clade rooted at this taxon','Number of reads assigned directly to this taxon','Rank code','NCBI taxonomy ID','Scientific name']
        csv_table.drop ([1])
        csv_table.set_index('Scientific name', inplace=True) #remove index colums
#convert tsv to csv first        
        csv_table.to_csv('kraken.csv.tmp' )
        
        #print  (csv_table)
        table = ""
        with open('kraken.csv.tmp', 'rb') as csvfile:
                
                reader = csv.reader(csvfile , delimiter=',')
                    #df.reset_index(drop=True, inplace=True)
                    #print (reader)

                table += "<table border='1'>\n"
                for row in reader:
                        table += "<tr>\n" + "".join(["<td>%s</td>\n" % 
                                     item for item in row]) + "</tr>\n"
                table += "</table>\n"
                #print(table)
                #print type (table)
	        # convert csv to html
 		with open("kraken.html", 'w') as htmlfile:
                	 htmlfile.write(table)

#======================================================================================================================
with open(outputfileconcat , "w") as f :
	message1 = """<html>
<head>
<style>
    header{ padding: 10px;
        text-align: "centre";
            background-image: linear-gradient(midnightblue, white);
            height: 250px ; width=1600px;
            color: white;
            font-size: 20px;  font-family: "Arial", sans-serif ; font-weight: bold;}
        body { color : white, font-family: cursive; sans-serif; overflow: visible; margin-bottom: 200%;}
    h1   {color: white; align="centre"}
    h2   {color: midnightblue; align="left" ; font-family: "Arial", sans-serif }
    h3   {color: midnightblue; align="centre" ; font-family: "Arial", sans-serif}
    p    {color: midnightblue;  font-family: "Arial", sans-serif}
    
</style>
</head>



<header>
<h1> General report </h1>
        <h2> Outputfile name : </h2> """
	message2 = outputfileconcat
	
	message3 = """
<image src= "https://pbs.twimg.com/media/DuIvwQxW4AAcwie.png" width="300" height="80" alt="" title="" style="float: top right;margin:-80 800; >
</header>

<body>
    <section id = "mapstat">
    <h2> Alignement results againt viral genomes </h2>
    <section id="hits">
    <embed src='Table.html'/  style="overflow:hidden;
    display:block; position:initial; height: 100%; width: 100%" >
    </section>

    <section id="kraken">
    <h2> kraken report  </h2>
    <embed src='kraken.html' style="overflow:hidden;
    display:block; position: initial; height: 100%; width: 100%" allowfullscreen="true">
    </section>

    <section id="krona">
    <h2> Interactive vizualisation of taxa relative abundance by Krona </h2>
    <iframe src="krona.html" frameborder="0"
    style=" display:inline; position: absolute; height: 100%; width: 100%"> </iframe>
    </section>

</body>
</html>

"""
	f.write(message1 + message2+ message3)
	f.close()

filenname =  outputfileconcat

webbrowser.open_new_tab(filenname)

