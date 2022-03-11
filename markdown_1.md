------------------------------------------------- pre orthology inference

Generate annotationÃ:
```
for i in $(ls raw_cds_reduced/*nt.fa); 
	do 
	base=$(echo $i | awk -F "/" '{print $NF}' | awk -F "_" '{print $1}'); 
	name=$(echo $i | awk -F "/" '{print $NF}' | awk -F "_" '{print $2}'); 
	echo -e "$base\t$name"; 
done > annotation.txt
```

Include species name 
```
for i in *fa; 
	do 
	b=$(echo $i | awk -F ".fa" '{print $1}'); 
	sed -e "s/^>/>$b|/" $i > $i.ref; 
done
```

Translate CDS to amminoacids
```
for i in *nt.fa; 
	do 
	base=$(echo $i | awk -F ".nt.fa" '{print $1}'); 
	transeq -table 11 -trim -sequence $i -outseq $base.aa.fa; 
done
```

------------------------------------------------- pre phylogenetic inference

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

Rretrotranslate OGs
```
for i in $(cat single_over168sp_og.lst); 
	do 
	for j in $(grep ">" $i); 
		do 
		sp=$(echo $j | awk -F "%" '{print $1}'| tr -d ">"); 
		gn=$(echo $j | awk -F "%" '{$1=""}1' | tr -d ">" | sed "s/_1$//" | sed "s/lcl /lcl%/"); 
		grep -w -A 1 $gn ../../../../../raw_cds/$sp*; 
	done > ../Orthogroup_Sequences_nt/$i.nt; 
done
```

MACSE preliminary alignment:
```
for i in *fa; do macse -prog alignSequences -seq $i --threads 16 -gc_def 11; done
```

Remove genes with inframe non-terminal stopcodons:
```
for i in *AA*; 
	do 
	for j in $(egrep -A 1 -e "\*-{,30}[A-Z]" -e "\?-{,30}[A-Z]" $i | grep ">"); 
		do 
		sed -i -e "/$j/,+1d" ../2_single_copy_min80%_nostop/${i::-8}NT.fasta; 
	done; 
done
```

Remove genes with frameshifts:
```
for i in *NT*; 
	do 
	awk 'NR==FNR{if (/!/) del[NR-1]; next} !(FNR in del)' $i $i | grep -v "\!"; 
done
```

Gblock nucleotides:
```
for i in *n.aln; do Gblocks $i -t=c -b5=h; done
```

Gblock amminoacids:
```
for i in *p.aln; do Gblocks $i -t=p -b5=h; done
```

Reformat after Gblocks:
```
for i in *gb; 
	do awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < $i > ${i::-6}gb.aln;
done
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






Sanity check:
```
for k in *NT.fasta; do echo -e "\n $k \n"; for i in $(grep ">" $k); do sp=$(echo $i | tr -d ">" | sed 's/\.1//' | sed 's/\.2//'); seq=$(grep $i -A 1 $k | tail -1 | tr -d "-"); hit=$(grep -w $seq ../../tot_formatted_nt/* | awk -F ".nt.fa" '{print $1}' | awk -F "/" '{print $NF}'); if echo $hit | grep -q $sp; then echo ok; else echo $sp $hit; fi; done; done
```

