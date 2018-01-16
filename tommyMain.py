from Bio import Phylo
from Bio.Phylo import PhyloXML
from Bio.Phylo.Applications import PhymlCommandline
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor, ParsimonyScorer, NNITreeSearcher, \
	ParsimonyTreeConstructor
from Bio.Align.Applications import MuscleCommandline
from Bio import AlignIO
from Bio import SeqIO

lookup = {}
with open("seqs.fasta", "rU") as handle:
    for record in SeqIO.parse(handle, "fasta"):
        lookup[record.id] = str(record.seq)
		

cline = MuscleCommandline(input="seqs.fasta", out="outAlign.aln", clw=True)
# cline()
# AlignIO.convert("outAlign.aln", "clustal", "outTree.phy", "phylip-relaxed")
#
# cmdline = PhymlCommandline(input='outTree.phy', datatype='aa', model='WAG', alpha='e', bootstrap=20)
# out_log, err_log = cmdline()
#
# egfr_tree = Phylo.read("outTree.phy_phyml_tree", "newick")
# Phylo.draw_ascii(egfr_tree)
#
# egfr_phy = egfr_tree.as_phyloxml()


#
stdout, stderr = cline()
# print(stdout, stderr)
# print(cline)
aln = AlignIO.read('outAlign.aln', 'clustal')
# print(aln)

calculator = DistanceCalculator('identity')
dm = calculator.get_distance(aln)
print(dm)

constructor = DistanceTreeConstructor(calculator,'nj')
tree = constructor.build_tree(aln)
print(tree)

scorer = ParsimonyScorer()
searcher = NNITreeSearcher(scorer)
constructor = ParsimonyTreeConstructor(searcher, tree)
pars_tree = constructor.build_tree(aln)

egfr_phy = pars_tree.as_phyloxml()

print(pars_tree)
print(egfr_phy)
print(list(pars_tree.find_elements('Inner3')))

for clade in egfr_phy.get_terminals():
    key = clade.name
    accession = PhyloXML.Accession(key, 'NCBI')
    mol_seq = PhyloXML.MolSeq(lookup[key], is_aligned=True)
    sequence = PhyloXML.Sequence(type='aa', accession=accession, mol_seq=mol_seq)
    clade.sequences.append(sequence)

Phylo.write(egfr_phy, 'egfr-family-annotated.xml', 'phyloxml')
