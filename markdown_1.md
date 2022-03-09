Find single-copy genes orthogroups:
```
for i in *fa; 
	do 
	a=$(grep ">" $i | wc -l); 
	b=$(grep ">" $i | awk -F "%" '{print $1}' | sort -u | wc -l); 
	if [[ $a == $b ]]; 
		then echo $i; 
	fi; 
done > single_og.lst
```

Find genes present in more than 80% of the species:
```
for i in $(cat single_og.lst); 
	do 
	n=$(grep -c ">" $i); 
	if [[ $n -gt 190 ]]; 
		then echo $i; 
	fi; 
done > single_over80%sp_og.lst
```

MACSE preliminary alignment:
```
for i in *fa; do macse -prog alignSequences -seq $i --threads 16 -gc_def 11; done
```

Remove genes with inframe non-terminal stopcodons:
```
for i in *AA*; do for j in $(egrep -A 1 "\*-{,30}[A-Z]" $i | grep ">"); do sed -i "/$j/{n;d}" ../2_single_copy_min80%_nostop/${i::-8}NT.fasta; done; done
```

Remove genes with frameshifts:
```
for i in *NT*; do awk 'NR==FNR{if (/!/) del[NR-1]; next} !(FNR in del)' $i $i | grep -v "\!"; done
```

Reformat after Gblocks:
```
for i in *gb; do awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < $i > ${i::-6}gb.aln; done
```

Concatenate the nucleotide alignements:
```
AMAS.py concat -i 
 -f fasta -d dna -t .fa  -y raxml -p .prt
```

Concatenate the amminoacid alignments:
```
AMAS.py concat -i
 -f fasta -d dna -t .fa  -y raxml -p .prt
```

RY recoding for nucleotide alignemnts:
```
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' single_copy_min80%_gb_n.fa | 
sed "1~2 s/A/R/g" | sed '1~2 s/G/R/g' | 
sed "1~2 s/C/Y/g" | sed '1~2 s/T/Y/g'
```

Dayhoff-6 recoding for amminocid alignements:
```
```

