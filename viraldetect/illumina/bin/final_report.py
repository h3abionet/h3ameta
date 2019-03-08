import pandas as pd
import webbrowser
import os
import sys 

def main(argv):
    if len(argv) >= 2:
        print ("\nUsage: python2.7 final_report.py <krakenfile> <mappingstatoutput> <kronaoutput> <sample> <outputfile>")
        sys.exit(2)

krak = sys.argv[1]
stats =  sys.argv[2]
krona =  sys.argv[3]
sample =  sys.argv[4]+" sample"
outputfileconcat = sys.argv[5]

def csv_to_html_table(fname,headers=None,delimiter=",", tag=""):
    with open(fname) as f:
        content = f.readlines()
    #reading file content into list
    rows = [x.strip() for x in content] 
    table = "<table class=\"scroll-table\""+  "height=200" + " id=\"" + tag + "\">"
    #print(table)
    #creating HTML header row if header is provided 
    if headers is not None:
        table+= "".join(["<th>"+cell+"</th>" for cell in headers.split(",")])
    else:
        table+= "".join(["<th>"+cell+"</th>" for cell in rows[0].split(delimiter)])
        rows=rows[1:]
    #Converting csv to html row by row
    for row in rows:
        table+= "<tr>" + "".join(["<td>"+cell+"</td>" for cell in row.split(delimiter)]) + "</tr>" + "\n"
    table+="</table><br>"
    return table

col_head_map='Pathogen, Genome size, Number of mapped reads, Number of Unmapped reads'
mapStats=csv_to_html_table(stats, delimiter="\t", headers=col_head_map, tag="t02")
col_headers='Reads rooted at taxon(%), Reads rooted at this taxon,Reads assigned to taxon,Rank code, Taxon ID, Scientific Name'
krak_res = csv_to_html_table(krak, delimiter="\t", headers=col_headers, tag="t01")

with open(outputfileconcat , "w") as f :
    message1 = '<html><head><header><h1 align="center"> Viral Detection report for '+sample+'</h1>'
    message3 = """
    <head>
    <title> Input file name </title>
    <script src="https://ajax.googleapis.com/ajax/libs/jquery/3.1.0/jquery.min.js"></script>
    <script type="text/javascript" src="https://code.jquery.com/jquery-3.3.1.js"></script>
    <script  type="text/javascript" src="https://cdn.datatables.net/1.10.19/js/jquery.dataTables.min.js"></script>
    <script type="text/javascript" src="https://cdnjs.cloudflare.com/ajax/libs/jszip/3.1.3/jszip.min.js"></script>
    <script type="text/javascript" src="report.js"></script>

    <link rel="stylesheet" type="text/css" href="https://cdn.datatables.net/1.10.19/css/jquery.dataTables.min.css">
    <style type="text/css">
    #t01 {
      font-family: "Trebuchet MS", Arial, Helvetica, sans-serif;
      border-collapse: collapse;
      width: 80%;
    }

    #t01 td, #t01 th {
      border: 1px solid #ddd;
      padding: 8px;
    }

    #t01 tr:nth-child(even){background-color: #D8E0E5;}

    #t01 tr:hover {background-color: #ddd;}

    #t02 th {
      padding-top: 12px;
      padding-bottom: 12px;
      text-align: left;
      background-color: #417AC6;
      color: white;
    }
    #t02 {
      font-family: "Trebuchet MS", Arial, Helvetica, sans-serif;
      border-collapse: collapse;
      width: 80%;
    }

    #t02 td, #t01 th {
      border: 1px solid #ddd;
      padding: 8px;
    }

    #t02 tr:nth-child(even){background-color: #D8E0E5;}

    #t02 tr:hover {background-color: #ddd;}

    #t02 th {
      padding-top: 12px;
      padding-bottom: 12px;
      text-align: left;
      background-color: #417AC6;
      color: white;
    }
    #mapping {background-color: #f2f2f2;}
    .scroll-table-container {height: 500px; overflow: scroll;background-color: #f2f2f2;}
    .scroll-table, td, th{border-collapse:collapse; border:1px solid #777; max-width: 700px;}
    </style>
    </head>
    </header>"""
    
    body="""<body bgcolor="#D8E0E5">
    <div align="center">
    <section id="mapping">
    <h2> Mapping Statistics </h2>
    """

    krona_bit="""
    </section>
    <h2> Visualisation of kraken results </h2>
        <br>
        <section id="hits">
        <div style="height: 700px;">
            <embed src='"""+krona+"""' width="800px" height="700px" />
        <div>
        </section>
    <div class="scroll-table-container">
    <h2>Tabular kraken results </h2>
    """
    
    last="""</div></section></div>
    <div align="center" style="width=10px">
    <p><font color="red"><Strong>
    DISCLAIMER: This pipeline is still under development, we hope to add more details regarding the process of
    generating these results, interpretation of results and Visualisation. 
    </Strong></font>
    </p>
    </div>
    <script src="ddff.js"></script>
    <script>
        $('#t01').ddTableFilter();
    </script>
    </body></html>"""
    f.write(message1 + message3 + body + mapStats + krona_bit + krak_res  + last)
    f.close()
filenname =  outputfileconcat
webbrowser.open_new_tab(filenname)

