#!/usr/bin/env python2.7
import os
import sys
import re
import glob
import csv
import pandas
from pandas import DataFrame
import numpy as np
import webbrowser
import argparse
import shutil, os

outpath = "./output"

kraken = sys.argv[1]
stat =  sys.argv[2]
krona =  sys.argv[3]
outputfileconcat = sys.argv[4]
search_pattern_path = sys.argv[5]

#==================================================================================================================
#Detecttion of viruses
stringToMatch= []
#search_pattern.txt file genereated from grep ">" all.fasta & cut -d "|" -f4,5hader_all_fasta_gerrit.txt > search_pattern.txt
with open(os.path.join(search_pattern_path), "r") as f:
    for s in f:
        #print(s)
        splitElements = s.split('|')
        id = splitElements[0]
        sp = splitElements[4]
        #print(sp)
        stringToMatch.append(sp)
#print type(stringToMatch)

df = []
for key in stringToMatch :
    splitElements = key.split(' ')
    key2= splitElements[1]
    krakenfile = open (kraken , "r")
    trouv = key,"Not Found"
    for line in krakenfile :
        #print(line)
        if key2 in line:
            trouv = (key,"Found")
    #print(trouv)
    df.append(trouv)
#print(df)
kr_df = DataFrame.from_records(df)
#print (kr_df)


#===================================================================================
#bowtie2  output

stat_df = pandas.read_csv(stat, sep="\t", header=None)
stat_df.columns = ["a", "b", "c", "d"]
stat_df = stat_df[:-1] #delete last line
stat_df["a"]= stat_df.a.str.split('|').str[3] #keep only NC
#print (stat_df)


#merge results
result = pandas.concat([kr_df, stat_df], axis=1)
#print(result)

result.columns = [ 'Pathogen', 'kraken result', 'Reference sequence identifier' ,'Reference sequence length', 'Number of mapped reads',  'Number of unmapped reads']
here = os.path.dirname(os.path.realpath(__file__))
subdir = "../output"
filenametable = "Table.html"
filepathtable = os.path.join(here, subdir, filenametable)
result.to_html(filepathtable, index=False)


#=========================================================================================================================
#convert kraken tsv to html :#https://github.com/dbohdan/csv2html?fbclid=IwAR0L3a41xXXVpWXxH-q-httRuB1x-6S6KXQysBAx5dGhBaE3Ntg6
#add header
table = ""
nargs = len(sys.argv)
with open(kraken, 'rb') as tsvfile :
        csv_table=pandas.read_csv(tsvfile , sep='\t')
        csv_table.columns = ['Percentage of reads covered by the clade rooted at this taxon' ,'Number of reads covered by the clade rooted at this taxon','Number of reads assigned directly to this taxon','Rank code','NCBI taxonomy ID','Scientific name']
        csv_table.drop ([1])
        csv_table.set_index('Scientific name', inplace=True) #remove index colums
#convert tsv to csv first
	filenamekr = "kraken.csv.tmp"
	filepathkr = os.path.join(here, subdir, filenamekr)
        csv_table.to_csv(filepathkr)

        #print(csv_table)
        table = ""
        with open(filepathkr, 'rb') as csvfile:

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
		filenamekrHTML = "kraken.html"
        	filepathkrHTML = os.path.join(here, subdir, filenamekrHTML)
 		with open(filepathkrHTML, 'w') as htmlfile:
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
        <h2> Outputfile name : </h2>
 """
    message2 = outputfileconcat

    message3 = """

        <image src= "https://pbs.twimg.com/media/DuIvwQxW4AAcwie.png" width="300" height="80" alt="" title="" style="float: top right;margin:-80 800; >
        </header>

        <body>
        <section id = "mapstat">
        <h2> Alignement results againt viral genomes </h2>
        <section id="hits">
        <embed src= "./Table.html"  style="overflow:hidden;
        display:block; position:initial; height: 100%; width: 100%" >
        </section>

        <section id="kraken">
        <h2> kraken report  </h2>
        <embed src="./kraken.html" style="overflow:hidden;
        display:block; position: initial; height: 100%; width: 100%" allowfullscreen="true">
        </section>

        <section id="krona">
        <h2> Interactive vizualisation of taxa relative abundance by Krona </h2>
        <iframe src="./krona.html" frameborder="0"
        style=" display:inline; position: absolute; height: 100%; width: 100%"> </iframe>
        </section>

        </body>
        </html>

        """
    f.write(message1 + message2+ message3)
    f.close()

filenname =  outputfileconcat
#webbrowser.open_new_tab(filenname)
