# projet_python

Python 3.7 recommandé pour l'éxecution des script

Dépendances : Pandas / Numpy / Matplotlib

Le projet comporte une librarie principale (main_library.py) ainsi que 3 scripts (split_chains.py, find_salt_bridges.py, contact_map.py)

## La librairie principale

Cette libraire est composée de 9 fonctions distinctes :

### read_pdb()

Read_pdb est une fonction qui à partir d'un fichier .pdb, créer un dataframe pandas comportant toutes les coordoonées des atomes et hétéroatomes définis par "ATOM" et "HETATM". Cette fonction crée une liste pour chaque coordonnée, ces listes sont passées à un dictionnaires ordonnées qui est à son tour transformé en dataframe pandas.

### write_pdb()

Write_pdb est une fonction qui à partir d'un dataframe pandas composé des coordonnées d'un fichier .pdb, écris un nouveau fichier avec ces coordonnées dans le format .pdb. Cette fonction prend en input le nom du fichier de sortie ainsi que le dataframe pandas.

### select_atoms()

Select_atoms est une fonction qui retourne un sous-ensemble d'un dataframe, à partir d'un dataframe pandas et de critères de sélections sous la forme d'un dictionnaire. Le dictionnaire doit avoir les valeurs associées aux clés sous la forme de liste. Un dataframe temporaire est créer pour éviter de modifier le dataframe de base. Pour chaque clés du dictionnaire, toutes les lignes du dataframe de base correspondant à ses valeurs sont jointes à un dataframe vide.

### split_chains()

Split_chains est une fonction qui crée un fichier de coordonnées des atomes par chaîne, dans le format .pdb. Cette fonction prend comme input un fichier dans le format .pdb. Les nouveaux fichier .pdb produit sont sous la forme 'nompdb_nomchaîne'. Un .pdb est crée par chaînes distinctes dans le .pdb initial.

### aa321()

aa321 est une fonction qui à partir d'une liste d'acides aminés à 3 lettres, retournes une liste de séquences d'acides aminés à 1 lettres. 

### get_aa_seq()

get_aa_seq est une fonction qui prends en argument un dataframe pandas, et retourne une liste composé de la séquence d'acide aminés à 1 lettres. Chaque chaine est séparé en chaines de caractère distinctes. Cette fonction extrait la séquence d'acides aminés 3lettres pour chaques chaines et la convertie en séquence 1 lettres par la fonction split_chains().

### compute_distance()

compute_distance est une fonction qui prend deux lignes d'un dataframe pandas, composés de coordonnées atomique d'un pdb, et calcule la distance euclidienne entre ces deux lignes representant deux atomes. Cette fonction calcul la distance euclidienne à partir des coordonnées x,y et z retrouvés dans chaque lignes. Si une de ces distances n'est pas présentes dans la ligne, la fonction renvoie une erreur.

### find_salt_bridges()

Find_salt_bridges est une fonction qui prend comme argument un dataframe, composé des coordonnées atomiques d'un pdb ainsi qu'une valeur sous la forme d'un entier répresentant le cutoff. Cette fonction mesures toutes les distances deux à deux entres les atomes donneurs et accepteurs des acides aminés chargés. Si cette distance est en-dessous de la valeur du cutoff, la fonction considère la présence d'un pont salins et renvoie la paires d'acides aminés engagée dans cette liaison. Toutes les paires déterminés sont stockées dans une liste.

### contact_map()



