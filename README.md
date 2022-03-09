
```
for i in *fa; 
do a=$(grep ">" $i | wc -l); b=$(grep ">" $i | awk -F "%" '{print $1}' | sort -u | wc -l); 
if [[ $a == $b ]]; 
then echo $i; 
fi; 
done > single_og.lst
```

