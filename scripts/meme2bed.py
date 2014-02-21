from Bio import motifs
import csv

handle = open('meme.txt', 'r')
memerecord = motifs.parse(handle, 'meme')
handle.close()

for motif in memerecord:
	bed = []
	for feature in motif.instances:
		chromossome = feature.sequence_name.split(':')[0]
		start = int(feature.sequence_name.split(':')[1].split('-')[0]) + feature.start
		end = start + motif.length
		bed.append(
			[
				chromossome,
				start,
				end,
				motif.consensus,
				feature.pvalue,
				feature.strand
			]
		)
	file_name = motif.name.replace(" ", "_") + '.bed'
	resultFile = open(file_name,'wb')
	wr = csv.writer(resultFile, delimiter='\t', quoting=csv.QUOTE_MINIMAL)
	wr.writerows(bed)
	resultFile.close()