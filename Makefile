##############################
# Makefile du Projet
##############################

#Modifier le chargement des fichiers
OBJECTSMAIN = Element.o Atom.o PeriodicTable.o Structure.o Domain.o Grid.o main.o
all: main.exe

#Commande pour la compilation et le chargement de ROOT
COMPILCOMMAND=g++ -c -o $@   $<

# Commande generique pour generer tous les fichiers objets
# a partir des fichiers d'en tete et des fichiers source
%.o: %.cpp %.h
	$(COMPILCOMMAND)

# La regle ci-dessus ne marche pas necessairement pour main.o du fait de
# l'absence d'un fichier main.h...
%.o: %.cpp 
	$(COMPILCOMMAND)

# Fichiers source a generer
# A noter les "dependances" en plus des fichiers d'en tete et des fichiers
# source definis dans la commande generique ci-dessus
Element.o:		Element.h
Atom.o:			Atom.h		Constants.h
PeriodicTable.o: 	PeriodicTable.h
Structure.o:		Structure.h
Domain.o:		Domain.h
Grid.o:			Grid.h		Constants.h
main.o:

# Creation de l'executable
main.exe: $(OBJECTSMAIN)
	g++ -o main.exe $(OBJECTSMAIN)
	

# Commande de nettoyage
clean:
	rm -rf *.o main.exe *~ Print
