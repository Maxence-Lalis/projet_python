# projet_python

Python 3.7 recommandé pour l'éxecution des script

Dépendances :
Pandas
Numpy
Matplotlib

Le projet comporte une librarie principale (main_library.py) ainsi que 3 scripts (split_chains.py, find_salt_bridges.py, contact_map.py)

## La librairie principale

Cette libraire est composée de 10 fonctions distinctes :

### read_pdb()

Read_pdb est une fonction qui à partir d'un fichier .pdb, créer un dataframe pandas comportant toutes les coordoonées des atomes et hétéroatomes définis par "ATOM" et "HETATM". Cette fonction crée une liste pour chaque coordonnée, ces listes sont passées à un dictionnaires ordonnées qui est à son tour transformé en dataframe pandas.

### write_pdb()

Write_pdb est une fonction qui à partir d'un dataframe pandas composé des coordonnées d'un fichier .pdb, écris un nouveau fichier avec ces coordonnées dans le format .pdb. Cette fonction prend en input le nom du fichier de sortie ainsi que le dataframe pandas.

### select_atoms()

Select_atoms est une fonction qui retourne un sous-ensemble d'un dataframe, à partir d'un dataframe pandas et de critères de sélections sous la forme d'un dictionnaire. Le dictionnaire doit avoir les valeurs associées aux clés sous la forme de liste. Un dataframe temporaire est créer pour éviter de modifier le dataframe de base. Pour chaque clés du dictionnaire, toutes les lignes du dataframe de base correspondant à ses valeurs sont jointes à un dataframe vide.

### split_chains()

Split_chains est une fonction qui crée un fichier de coordonnées des atomes par chaîne, dans le format .pdb. Cette fonction prend comme input un fichier dans le format .pdb. Les nouveaux fichier .pdb produit sont sous la forme 'nompdb_nomchaîne'. Un .pdb est crée par chaînes distinctes dans le .pdb initial.

### aa321()

aa321 est une fonction qui à partir d'une liste d'acides aminés à 3 lettres, retournes une liste de séquences d'acides aminés à 1 lettres par chaînes. 

### get_aa_seq()

