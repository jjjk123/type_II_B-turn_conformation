# http://thegrantlab.org/bio3d/reference/torsion.pdb.html#arguments
install.packages("bio3d", dependencies=TRUE)

pdb <- read.pdb("1g60")
tor <- torsion.pdb(pdb)

write.csv(tor$tbl, '1g60_angles.csv')
