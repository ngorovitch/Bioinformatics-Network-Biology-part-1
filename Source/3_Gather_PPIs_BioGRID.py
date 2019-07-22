import pandas as pd
import csv

# Funtion that returns the list of interaction tuples
# Parameters:
#   PPI_data: PPI source filtered data
#   protein_list: list of seed proteins (e.g. if we intend to filter by 'GeneID' it should be a list of 'GeneID')
#   DB_protein_A_col_name_search: the name of the column that we intend to use to perfom the search
#   the output is in the following form: ([(p1_a, p2_b), (p3_a, p4_b), ...], [(p1_c, p2_d), (p3_c, p4_d), ...])
#   where a and c represent two different types of protein ID standing for the same protein.
#   DB_protein_A_col_name1: The name of the column for the first type of protein ID and for the first protein
#   DB_protein_B_col_name1: The name of the column for the first type of protein ID and for the second protein
#   DB_protein_A_col_name2: The name of the column for the second type of protein ID and for the first protein
#   DB_protein_B_col_name2: The name of the column for the second type of protein ID and for the second protein
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
    #### ------------Perform PPI search using BioGRID---------------
    
    ## ------Searching by "GeneID"--------
    
    print()
    print('Performing the search of PPIs using BioGRID source...')
    # Load the full BioGRID PPI source data
    BioGRID = pd.read_csv('../Data/BioGRID/BIOGRID-ORGANISM-Homo_sapiens-3.5.167.tab2.txt', sep = '\t',  index_col=0, dtype={'Score': str})    
    # Load our previously created and manually integrated dataset
    GenesInfo = pd.read_csv('../Question 2/Manually_integrated_basic_info.csv', sep = ',',  index_col=0)
    # To make things fast we are going to filter our data by the 'Entrez Gene Interactor A'
    # we want to remain only with the rows of our interest
    BioGRID_seed1 = BioGRID.loc[BioGRID['Entrez Gene Interactor A'].isin(GenesInfo['Entrez Gene ID'])]
    # This list (our interactome) is going to contain all the protein protein interactions as a list of tuples
    # [(protein_1_A, protein_1_B), (protein_2_A, protein_2_B), (protein_1_a, protein_1_b), ...]
    All_PPIs_bioGRID1 = get_interactome(BioGRID_seed1, GenesInfo['Entrez Gene ID'], 'Entrez Gene Interactor A', 
                                        'Entrez Gene Interactor A', 'Entrez Gene Interactor B',
                                        'Official Symbol Interactor A', 'Official Symbol Interactor B')
    # Let's retrieve the list of all the non-seed proteins interacting with at least one seed gene
    non_seed_Proteins1 = list(set([x for x in [e[1] for e in All_PPIs_bioGRID1] if x not in list(GenesInfo['Entrez Gene ID'])]))
    # Let's retrieve and include in our interactome the interactions among these non-seed proteins
    BioGRID_non_seed1 = BioGRID.loc[BioGRID['Entrez Gene Interactor A'].isin(non_seed_Proteins1)]
    BioGRID_non_seed1 = BioGRID.loc[BioGRID['Entrez Gene Interactor B'].isin(non_seed_Proteins1)]
    All_PPIs_bioGRID1 += get_interactome(BioGRID_non_seed1, non_seed_Proteins1, 'Entrez Gene Interactor A', 
                                        'Entrez Gene Interactor A', 'Entrez Gene Interactor B',
                                        'Official Symbol Interactor A', 'Official Symbol Interactor B')
    # remove duplicated PPIs
    All_PPIs_bioGRID1 = sorted(list(set(All_PPIs_bioGRID1)), key=lambda tup: tup[0])
    
    
    ## ------Searching by "gene symbols"--------
    
    # To make things fast we are going to filter our data by the 'Entrez Gene Interactor A'
    # we want to remain only with the rows of our interest
    BioGRID_seed2 = BioGRID.loc[BioGRID['Official Symbol Interactor A'].isin(GenesInfo['HGNC Approved gene symbol'])]
    # This list (our interactome) is going to contain all the protein protein interactions as a list of tuples
    # [(protein_1_A, protein_1_B), (protein_2_A, protein_2_B), (protein_1_a, protein_1_b), ...]
    All_PPIs_bioGRID2 = get_interactome(BioGRID_seed2, GenesInfo['HGNC Approved gene symbol'], 'Official Symbol Interactor A',
                                       'Entrez Gene Interactor A', 'Entrez Gene Interactor B',
                                        'Official Symbol Interactor A', 'Official Symbol Interactor B')
    # Let's retrieve the list of all the non-seed proteins interacting with at least one seed gene
    non_seed_Proteins2 = list(set([x for x in [e[3] for e in All_PPIs_bioGRID2] if x not in list(GenesInfo['HGNC Approved gene symbol'])]))
    # Let's retrieve and include in our interactome the interactions among these non-seed proteins
    BioGRID_non_seed2 = BioGRID.loc[BioGRID['Official Symbol Interactor A'].isin(non_seed_Proteins2)]
    BioGRID_non_seed2 = BioGRID.loc[BioGRID['Official Symbol Interactor B'].isin(non_seed_Proteins2)]
    All_PPIs_bioGRID2 += get_interactome(BioGRID_non_seed2, non_seed_Proteins2, 'Official Symbol Interactor A',
                                       'Entrez Gene Interactor A', 'Entrez Gene Interactor B',
                                        'Official Symbol Interactor A', 'Official Symbol Interactor B')
    # remove duplicated PPIs
    All_PPIs_bioGRID2 = sorted(list(set(All_PPIs_bioGRID2)), key=lambda tup: tup[0])
    
    # For each DB, check if the results are different when carrying out the search by using gene
    # symbols or by using Uniprot AC identifiers (in this case we used GeneID). Discrepancies must be reported.
    if(not All_PPIs_bioGRID1 == All_PPIs_bioGRID2):
        print()
        print('     Discrepancies observed between search by \'GeneIDs\' an search by \'gene symbols\':')
        print('         Search by GeneID: '+str(len(All_PPIs_bioGRID1))+' interactions')
        print('         Search by gene symbol: '+str(len(All_PPIs_bioGRID2))+' interactions')
    
    # Store the data gathered into a csv file
    with open('../Question 3/PPI_BioGRID.csv','w') as out:
        csv_out=csv.writer(out, lineterminator='\n')
        csv_out.writerow(['', 'Entrez Gene Interactor A','Entrez Gene Interactor B', 'Official Symbol Interactor A', 'Official Symbol Interactor B'])
        for tup in All_PPIs_bioGRID2:
            csv_out.writerow(tuple([0]+list(tup)))
    
    # Summarize the results to put in a table
    # no. of seed genes found
    all_genes = list(set([tup[0] for tup in All_PPIs_bioGRID2] + [tup[1] for tup in All_PPIs_bioGRID2]))
    seed_genes = list(set([gene for gene in all_genes if gene in list(GenesInfo['Entrez Gene ID'])]))
    number_seed_genes_found = len(seed_genes)
    # total no. of interacting proteins, including seed genes
    total_number_of_proteins = len(all_genes)
    # total no. of interactions found
    total_number_of_interactions = len(All_PPIs_bioGRID2)
    
    print()
    print('     Summarize of the results:')
    print('         Number of seed genes found: ', number_seed_genes_found)
    print('         Total Number of interacting proteins: ', total_number_of_proteins)
    print('         Total Number of interactions found: ', total_number_of_interactions)
    