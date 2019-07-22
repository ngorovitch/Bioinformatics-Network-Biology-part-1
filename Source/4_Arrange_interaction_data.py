import pandas as pd
import csv

# We want to build and store 3 tables:
#    seed_genes_interactome.csv
#    union_interactome.csv
#    intersection_interactome.csv

if __name__ == "__main__":    
    #read the GeneID - UniprotAC mapping from the disc
    with open('../Data/GeneID_to_UniprotAC.csv', mode='r') as infile:
        reader = csv.reader(infile)
        GeneID_to_UniprotAC = {rows[0]:rows[1] for rows in reader}
    
    # Load our previously created and manually integrated dataset
    GenesInfo = pd.read_csv('../Question 2/Manually_integrated_basic_info.csv', sep = ',',  index_col=0)
    
    interactome_BioGRID = pd.DataFrame()
    # read sotored BioGRID PPI data
    PPI_BioGRID = pd.read_csv('../Question 3/PPI_BioGRID.csv', sep = ',',  index_col=0)
    interactome_BioGRID['interactor A gene symbol'] = PPI_BioGRID['Official Symbol Interactor A']
    interactome_BioGRID['interactor B gene symbol'] = PPI_BioGRID['Official Symbol Interactor B']
    interactome_BioGRID['interactor A Uniprot AC'] = PPI_BioGRID['Entrez Gene Interactor A'].apply(lambda GeneID: GeneID_to_UniprotAC[str(GeneID)].split(';')[0].strip())
    interactome_BioGRID['interactor B Uniprot AC'] = PPI_BioGRID['Entrez Gene Interactor B'].apply(lambda GeneID: GeneID_to_UniprotAC[str(GeneID)].split(';')[0].strip())
    interactome_BioGRID['database source'] = 'BioGRID'
    
    interactome_IID = pd.DataFrame()
    # read sotored IID PPI data
    PPI_IID = pd.read_csv('../Question 3/PPI_IID.csv', sep = ',',  index_col=0)
    interactome_IID['interactor A gene symbol'] = PPI_IID['Official Symbol Interactor A']
    interactome_IID['interactor B gene symbol'] = PPI_IID['Official Symbol Interactor B']
    interactome_IID['interactor A Uniprot AC'] = PPI_IID['Uniprot AC ID Interactor A']
    interactome_IID['interactor B Uniprot AC'] = PPI_IID['Uniprot AC ID Interactor B']
    interactome_IID['database source'] = 'IID'

    interactome = interactome_BioGRID.append(interactome_IID)
    
    #seed genes interactome: interactions that involve seed genes only, from all DBs, in
    #the format:
    #    interactor A gene symbol, interactor B gene symbol, interactor A Uniprot AC, interactor B
    #    Uniprot AC, database source
    seed_genes_interactome = interactome.loc[ (interactome['interactor A Uniprot AC'].isin(GenesInfo['Uniprot AC'])) &  
                                             (interactome['interactor B Uniprot AC'].isin(GenesInfo['Uniprot AC']))]
    seed_genes_interactome.to_csv('../Question 4/seed_genes_interactome.csv', sep = ',', encoding='utf-8', index = False)
    
    # union interactome: all proteins interacting with at least one seed gene, from all DBs,
    # same format as above.
    union_interactome = interactome.loc[ (interactome['interactor A Uniprot AC'].isin(GenesInfo['Uniprot AC'])) |  
                                             (interactome['interactor B Uniprot AC'].isin(GenesInfo['Uniprot AC']))]
    union_interactome.to_csv('../Question 4/union_interactome.csv', sep = ',', encoding='utf-8', index = False)
    
    #intersection interactome: all proteins interacting with at least one seed gene confirmed
    #by both DBs, in the format:
    #    interactor A gene symbol, interactor B gene symbol, interactor A Uniprot AC, interactor B
    #    Uniprot AC
    union_interactome_BioGRID = interactome_BioGRID.loc[ (interactome_BioGRID['interactor A Uniprot AC'].isin(GenesInfo['Uniprot AC'])) |  
                                             (interactome_BioGRID['interactor B Uniprot AC'].isin(GenesInfo['Uniprot AC']))]
    union_interactome_IID = interactome_IID.loc[ (interactome_IID['interactor A Uniprot AC'].isin(GenesInfo['Uniprot AC'])) |  
                                             (interactome_IID['interactor B Uniprot AC'].isin(GenesInfo['Uniprot AC']))]
    
    intersection_interactome = union_interactome_BioGRID.drop(columns = ['database source']).merge(union_interactome_IID.drop(columns = 'database source'))
    intersection_interactome.to_csv('../Question 4/intersection_interactome.csv', sep = ',', encoding='utf-8', index = False)