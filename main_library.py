# -*- coding: utf-8 -*-
"""
Created on Mon Oct 19 18:32:47 2020

@author: Max
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Fonction read_pdb : Retourne un dataframe pandas, en prenant comme input un fichier .pdb.
def read_pdb(filename):
   
    # Création des listes vides pour chaques coordonnées
    coord_x = []; coord_y = []; coord_z = []; res_number = []; res_name = []
    atom_name = [] ; chain_id = []; type_atom = []; alt_loc_ind = []
    atom_number = []; CIR = []; occupancy = []; T_factor = []; E_symbol = [] 
    charge = []
    
    # Pour chaque ligne commencant par "ATOM" ou "HETATM", les coordonnées sont stockées
    # dans les listes correspondantes
    with open(filename ,"r") as pdb:
        for line in pdb:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                type_atom.append(str(line[0:6].strip()))
                atom_number.append(int(line[6:11]))
                alt_loc_ind.append(str(line[16:17].strip()))
                CIR.append(str(line[26:27].strip()))
                coord_x.append(float(line[30:38]))
                coord_y.append(float(line[38:46]))
                coord_z.append(float(line[46:54]))
                res_number.append(int(line[22:26]))
                res_name.append(str(line[17:20].strip()))
                atom_name.append(str(line[12:16].strip()))
                chain_id.append(str(line[21].strip()))
                occupancy.append(float(line[54:60]))
                T_factor.append(float(line[60:66]))
                E_symbol.append(str(line[76:78].strip()))
                charge.append(str(line[78:80].strip()))
        
        # Création d'un dictionnaire comportant chacune des listes de coordonnées
        extract = {"Type":type_atom, "Atom_Serial_Number":atom_number, 
                   "Atom_Name":atom_name, 
                   "alternate_location_indicator":alt_loc_ind, 
                   "Residue_Name":res_name, "chain_id":chain_id,  
                   "Residue_Number":res_number, "Code_Insertion_Residues":CIR,
                   "X_Coordinates": coord_x, "Y_Coordinates":coord_y , 
                   "Z_Coordinates":coord_z, "Occupancy":occupancy,
                   "Temperature_Factor":T_factor, "Element_Symbol":E_symbol, 
                   "Charge":charge}
        
        # Transformation du dictionnaire en pandas DataFrame
        extract_dataframe = pd.DataFrame(extract)
        return(extract_dataframe)
    

# Fonction write_pdb: Ecris un fichier dans le format .pdb, en prenant comme
# input un dataframe
def write_pdb(new_file, dataframe):
    
    # Pour chaque ligne du dataframe, les coordonnées sont écrites dans le format
    # standard .pdb
    with open(new_file, "w+") as new:
        for i in dataframe.values:
            new.write("{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}{:1s}   "
                      "{:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          "
                      "{:2>s}{:2s}\n".format(*i))
        new.write("TER")
 

# Fonction select_atoms: Retourne un sous-ensemble de dataframe, à partir d'un dataframe
# de base et de critères de sélections sous la forme d'un dictionnaire
def select_atoms(df, selector):
    
    # Les valeurs du dictionnaires doivent être en format liste
    # Un dataframe temporaire est créer pour éviter de modifier le dataframe de base
    temp = df
    selection = pd.DataFrame()
    
    # Pour chaque clés du dictionnaire, toutes les lignes correspondant à ses valeurs
    # sont jointes à un dataframe vide
    for keys, values in selector.items():
        for i in range(len(values)):
            selection = selection.append(temp[temp[keys] == values[i]])
        temp = selection
        selection = pd.DataFrame()
    return(pd.DataFrame(temp))


 
def split_chains(pdbfile):
     
     with open(pdbfile ,"r") as pdb:
         pdbfile = pdbfile.split(".")
         
         select = [ line for line in pdb if (line.startswith("ATOM") or \
                                             line.startswith("HETATM"))]
         chains = { line[21] for line in select if line[21] in 
                   'ABCDEFGHIJKLMNOPQRSTUVWXY'}
         
         for chain in chains:
            out = pdbfile[0] + "_" + chain + ".pdb"
            with open(out, "w") as o:
                o.writelines([line for line in select if line[21] == chain])


def get_aa_seq(dataframe):
    c = [] ; d = [] ; f = [] ; e =[]
    
    chains = {chain for chain in dataframe["chain_id"]}
    resid = {resid for resid in dataframe["Residue_Number"]}
    for chain in chains:
        select_chain = dataframe[dataframe["chain_id"] == chain]
        for i in resid:
            select_resid = select_chain[select_chain["Residue_Number"] == i]
            if select_resid.empty or (select_resid.iloc[0,1] == "HETATM"):
                next
            else:
                c.append(select_resid.iloc[0,4])
                d = "".join(c)
                e = aa321(d)
        f.append(e)
        c = []
    return(f)
    

def aa321(d):
    aa = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
         'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
         'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
         'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}
    
    upper_seq= d.upper()
    single_seq=''
    for i in range(int(len(upper_seq)/3)):
        single_seq += aa[upper_seq[3*i:3*i+3]]
    return(single_seq)

def compute_distance(dataframe, row1,row2):
    
    df1 = dataframe[["X_Coordinates","Y_Coordinates","Z_Coordinates"]]
    
    if  df1.iloc[row1].all() and df1.iloc[row2].all() == True:
        
        a = np.array((dataframe["X_Coordinates"][row1],
                      dataframe["Y_Coordinates"][row1],
                      dataframe["Z_Coordinates"][row1]))
        b = np.array((dataframe["X_Coordinates"][row2],
                      dataframe["Y_Coordinates"][row2],
                      dataframe["Z_Coordinates"][row2]))
        
        dist = np.sqrt(np.sum((a-b)**2))
        return(dist)
    else:
        return(False)

def compute_distance_mean_resid(dataframe, residue_one, residue_two):
    
    dataframe_mean = dataframe.groupby("Residue_Number").mean()
    mean_coord = dataframe_mean[["X_Coordinates","Y_Coordinates",
                                 "Z_Coordinates"]]
    
    a = np.array((mean_coord["X_Coordinates"][residue_one],
                  mean_coord["Y_Coordinates"][residue_one],
                  mean_coord["Z_Coordinates"][residue_one]))
    b = np.array((mean_coord["X_Coordinates"][residue_two],
                  mean_coord["Y_Coordinates"][residue_two],
                  mean_coord["Z_Coordinates"][residue_two])) 
    
    dist_mean = np.sqrt(np.sum((a-b)**2))
    return(dist_mean)

def find_salt_bridges(dataframe, cutoff):
    
    pos = select_atoms(dataframe, {"Residue_Name":["ARG","LYS","HIS"], 
                                   "Atom_Name":["OE1","OE2","NH1","NH2","NZ",
                                                "NE2","OD1","OD2"]})
    neg = select_atoms(dataframe, {"Residue_Name":["ASP","GLU"]})
    
    salin = []
    for i in pos.values:
        a = np.array(i[8:11])
        for j in neg.values:
            b = np.array(j[8:11])
            dist = np.linalg.norm(a-b)
            posipair = "-".join([i[4],str(i[6])])
            negpair = "-".join([j[4],str(j[6])])
            if dist < cutoff and [negpair,posipair] not in salin:
                salin.append([negpair,posipair])
            else:
                next
    return(salin)

def contact_map(dataframe, choix, out):
    
    CA = select_atoms(dataframe, {"Atom_Name":["CA"],"chain_id":["A"]})
    ca_coord = CA[["X_Coordinates","Y_Coordinates","Z_Coordinates"]]
    
    matrix=np.zeros(((ca_coord.values).shape[0],(ca_coord.values).shape[0]))
    
    n_df=(ca_coord.values)
    
    for i in range((ca_coord.values).shape[0]):
        for j in range((ca_coord.values).shape[0]):
            matrix[i,j]=np.sqrt(np.sum((n_df[i]-n_df[j])**2))
    
    if choix:
        fig = plt.matshow(matrix, cmap = "rainbow")
        ax1 = fig.axes
        ax1.set_xlabel("Residue Number")
        ax1.xaxis.set_label_position('top')
        ax1.xaxis.tick_top()
        ax1.set_ylabel("Residue Number")
        plt.colorbar(label="distance")
        plt.savefig(out)
    else:
        matrix = matrix < 10
        fig = plt.matshow(matrix, cmap = "binary")
        plt.savefig(out)
