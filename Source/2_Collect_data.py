from bioservices import KEGG
from bioservices import UniProt
import pandas as pd
import re

'''
    This script collect and store the following basic information from a GeneList:
        • official gene symbol (check if the symbols are updated and approved on the HGNC
          website; report any issue/lack of data/potential misinterpretation)
        • Uniprot AC, ‘accession number’ (a.k.a. ’Uniprot entry’)
        • protein name (the main one only, do not report the aliases)
        • Entrez Gene ID (a.k.a. ‘GeneID’)
        • very brief description of its function (keep it very short, i.e. max 20 words)
        • notes related to the above information, if any and if relevant

'''

if __name__ == "__main__":
    # Load up the gene List
    geneList = "../Data/genes.txt"
    with open(geneList) as f:
        genes = [line.rstrip('\n') for line in f]
    
    # Initialization of the libraries to be used
    u = UniProt(verbose=False)
    k = KEGG(verbose=False)
    
    # Columns that will be kept after collection of data
    COLLUMNS = ['Gene names  (primary )', 'Entry', 'Protein names', 'Function [CC]','Status']
    
    # Final is our final dataframe. It aggregates append rows after each iteration 
    final= pd.DataFrame()
    for gene in genes:  
        try:    
            gene_name = gene
            print("Processing: ", gene)
            
            # Query the uniprot db to get information of interests. Sorted by score.
            # We are interested in Human Organisms
            query = 'gene:'+gene+' AND organism:"Homo sapiens (Human) [9606]"'
            res = u.search(query, frmt = "tab", columns = "id")
            
            # Fetch the most relevant result
            entry = res.split("\n")[1].split("\t")[0]
            print("Processing Uniprot ID: ", entry)
            # Access the details of the most relevant result
            df = u.get_df(re.sub('Entry', '', entry),1, "Homo sapiens (Human)")
            
            # Keep the columns of interests as cited in the Project description
            # This is summarized on top of the script
            df = df[COLLUMNS]
            
            # Get the GENE ID 
            gene_n = "up:" +entry
            res = k.conv("hsa",gene_n)
            if(type(res) == dict):
                df["Entrez gene id"] = re.sub('hsa:', '', res[gene_n])
                
            # If Gene ID no available, we will fetch replace by NA and fill it manually
            else:df["Entrez gene id"] = "NA"
            
            # cleaning the various strings to a user friendly format
            df["Protein names"] = re.sub(r"\(([^()]+)\)", "", re.sub(r"\(([^()]+)\)", "", df["Protein names"][0]))
            df["Protein names"] = re.sub(r"\[([^\[\]]+)\]", "", re.sub(r"\[([^\[\]]+)\]", "", df["Protein names"][0]))
            df["Function [CC]"] = re.sub(r"FUNCTION: ", "", df["Function [CC]"][0]).split(".")[0] + '.'
            df["Function [CC]"] = re.sub(r"\(([^()]+)\)", "", re.sub(r"\(([^()]+)\)", "", df["Function [CC]"][0]))
            df["Function [CC]"] = re.sub(r"\[([^\[\]]+)\]", "", re.sub(r"\[([^\[\]]+)\]", "", df["Function [CC]"][0]))
            
            # Removing when necessary the suffix "_HUMAN" from the Entry name to get the official gene symbol
            #df["Entry name"] = df["Entry name"][0].split('_')[0]
            df["Gene"] = gene
            # appending row to the final dataset
            final = final.append(df.iloc[0])
        except:
            pass
    
    cols = ['Gene', 'Gene names  (primary )', 'Entry', 'Protein names', 'Entrez gene id', 'Function [CC]', 'Status']
    new_cols = ['Input gene name', 'official gene symbol', 'Uniprot AC', 'protein name', 'Entrez Gene ID', 'Function', 'Status']
    final = final [cols]
    final.columns = new_cols
    # Store basic infos in TSV file
    final.to_csv('../Question 2/basic_info.csv', sep = ',', encoding='utf-8', index = True)