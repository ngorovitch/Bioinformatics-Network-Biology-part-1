import pandas as pd
import csv

def get_interactome(PPI_data, protein_list, DB_protein_A_col_name_search, 
                              DB_protein_A_col_name1, DB_protein_B_col_name1, 
                              DB_protein_A_col_name2, DB_protein_B_col_name2):
    PPIs = []
    for GeneID in set(protein_list):
        GeneID_df = PPI_data.loc[PPI_data[DB_protein_A_col_name_search] == GeneID]
        for index, row in GeneID_df.iterrows():
            PPIs.append((row[DB_protein_A_col_name1], row[DB_protein_B_col_name1], row[DB_protein_A_col_name2].split(';')[0].strip(), row[DB_protein_B_col_name2].split(';')[0].strip()))
    return PPIs

if __name__ == "__main__":

    ## -----------Perform search using IID (Integrated Interactions Database)-----------
    
    ## ------Searching by "Uniprot AC ID"--------
    
    print()
    print('Performing the search of PPIs using IID source...')

    All_PPIs_IID_full1 = []
    All_PPIs_IID_full2 = []
    
    # Load our previously created and manually integrated dataset
    GenesInfo = pd.read_csv('../Question 2/Manually_integrated_basic_info.csv', sep = ',')
    
    # Read the IID PPI source data by chunck
    for IID in  pd.read_csv('../Data/IID/human_annotated_PPIs.txt', sep = '\t', low_memory=False, 
                            chunksize=20000, usecols = ['uniprot1', 'uniprot2', 'symbol1', 'symbol2']):
        # To make things fast we are going to filter our data by the 'Entrez Gene Interactor A'
        # we want to remain only with the rows of our interest
        IID_seed1_column1 = IID.loc[IID['uniprot1'].isin(GenesInfo['Uniprot AC'])]
        IID_seed1_column2 = IID.loc[IID['uniprot2'].isin(GenesInfo['Uniprot AC'])]
        IID_seed1 = IID_seed1_column1.append(IID_seed1_column2)
        #This list (our interactome) is going to contain all the protein protein interactions as a list of tuples
        # [(protein_1_A, protein_1_B), (protein_2_A, protein_2_B), (protein_1_a, protein_1_b), ...]
        All_PPIs_IID1 = get_interactome(IID_seed1, GenesInfo['Uniprot AC'], 'uniprot1',
                                        'uniprot1', 'uniprot2',
                                        'symbol1', 'symbol2')
        # Let's retrieve the list of all the non-seed proteins interacting with at least one seed gene
        non_seed_Proteins1 = list(set([x for x in [e[1] for e in All_PPIs_IID1] if x not in list(GenesInfo['Uniprot AC'])]))
        # Let's retrieve and include in our interactome the interactions among these non-seed proteins
        IID_non_seed1 = IID.loc[IID['uniprot1'].isin(non_seed_Proteins1)]
        IID_non_seed1 = IID.loc[IID['uniprot1'].isin(non_seed_Proteins1)]
        All_PPIs_IID1 += get_interactome(IID_non_seed1, non_seed_Proteins1, 'uniprot1',
                                        'uniprot1', 'uniprot2',
                                        'symbol1', 'symbol2')
        #remove duplicated PPIs
        All_PPIs_IID1 = sorted(list(set(All_PPIs_IID1)), key=lambda tup: tup[0])
        
        ## ------Searching by "gene symbols"--------
        
        # To make things fast we are going to filter our data by the 'Entrez Gene Interactor A'
        # we want to remain only with the rows of our interest
        IID_seed2 = IID.loc[IID['symbol1'].isin(GenesInfo['HGNC Approved gene symbol'])]
        #This list (our interactome) is going to contain all the protein protein interactions as a list of tuples
        # [(protein_1_A, protein_1_B), (protein_2_A, protein_2_B), (protein_1_a, protein_1_b), ...]
        All_PPIs_IID2 = get_interactome(IID_seed1, GenesInfo['HGNC Approved gene symbol'], 'symbol1',
                                        'uniprot1', 'uniprot2',
                                        'symbol1', 'symbol2')
        # Let's retrieve the list of all the non-seed proteins interacting with at least one seed gene
        non_seed_Proteins2 = list(set([x for x in [e[3] for e in All_PPIs_IID2] if x not in list(GenesInfo['HGNC Approved gene symbol'])]))
        # Let's retrieve and include in our interactome the interactions among these non-seed proteins
        IID_non_seed2 = IID.loc[IID['symbol1'].isin(non_seed_Proteins2)]
        IID_non_seed2 = IID.loc[IID['symbol2'].isin(non_seed_Proteins2)]
        All_PPIs_IID2 += get_interactome(IID_non_seed1, non_seed_Proteins1, 'symbol1',
                                        'uniprot1', 'uniprot2',
                                        'symbol1', 'symbol2')
        #remove duplicated PPIs
        All_PPIs_IID2 = sorted(list(set(All_PPIs_IID2)), key=lambda tup: tup[0])
        
        All_PPIs_IID_full1 += All_PPIs_IID1
        All_PPIs_IID_full2 += All_PPIs_IID2
        
    
    All_PPIs_IID_full1 = sorted(list(set(All_PPIs_IID_full1)), key=lambda tup: tup[0])
    All_PPIs_IID_full2 = sorted(list(set(All_PPIs_IID_full2)), key=lambda tup: tup[0])
    
    #For each DB, check if the results are different when carrying out the search by using gene
    #symbols or by using Uniprot AC identifiers. Discrepancies must be reported.
    if(not All_PPIs_IID_full1 == All_PPIs_IID_full2):
        print()
        print('     Discrepancies observed between search by \'Uniprot AC IDs\' an search by \'gene symbols\':')
        print('         Search by Uniprot AC ID: '+str(len(All_PPIs_IID_full1))+' interactions')
        print('         Search by gene symbol: '+str(len(All_PPIs_IID_full2))+' interactions')
        
    #Store the data gathered into a csv file
    with open('../Question 3/PPI_IID.csv','w') as out:
        csv_out=csv.writer(out, lineterminator='\n')
        csv_out.writerow(['', 'Uniprot AC ID Interactor A','Uniprot AC ID Interactor B', 'Official Symbol Interactor A', 'Official Symbol Interactor B'])
        for tup in All_PPIs_IID_full1:
            csv_out.writerow(tuple([0]+list(tup)))
            
    # Summarize the results to put in a table
    # no. of seed genes found
    all_genes = list(set([tup[0] for tup in All_PPIs_IID_full1] + [tup[1] for tup in All_PPIs_IID_full1]))
    seed_genes = list(set([gene for gene in all_genes if gene in list(GenesInfo['Uniprot AC'])]))
    number_seed_genes_found = len(seed_genes)
    # total no. of interacting proteins, including seed genes
    total_number_of_proteins = len(all_genes)
    # total no. of interactions found
    total_number_of_interactions = len(All_PPIs_IID_full1)
    
    print()
    print('     Summarize of the results:')
    print('         Number of seed genes found: ', number_seed_genes_found)
    print('         Total Number of interacting proteins: ', total_number_of_proteins)
    print('         Total Number of interactions found: ', total_number_of_interactions)