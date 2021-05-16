echo "# Started" > list_of_commands.txt
exp_unique_num=$(awk 'NR==1 {next} {print $1}' ../../Merged_Data_5.csv | cut -d "," -f 4 | uniq | wc -l)
traits=64
seed=5
K=10
obj_best=0.03
Iter=10000000
for ((trait=1; trait<=traits; trait=trait*2)); do
	for ((i=1; i<=$exp_unique_num; i=i+10)); do 
		if [ $(($exp_unique_num - $i)) -lt 10 ]; then 
			echo "Almost done"; 
			last_j=$exp_unique_num
			unset exp_unique_index
		else 
			last_j=$(($i+$K-1))
		fi
		declare -a exp_unique_index
		for ((j=i; j<=last_j; j++)); do
			exp_index=$(awk 'BEGIN {FS=","} NR==1 {next} $4 == "'$j'" {print NR-1}' ../../Merged_Data_5.csv)
			arr=($exp_index)
			start=${arr[0]}
			end=${arr[-1]}
			current_exp=$(($(($j-$i))*3))
			exp_unique_index[$current_exp]=$j
			exp_unique_index[$(($current_exp+1))]=$start
			exp_unique_index[$(($current_exp+2))]=$end
		done
		echo "cd ./seed_${seed}traits_$trait; Rscript Run_shell_C.R ../../../Merged_Data_5.csv $seed $trait $obj_best $Iter ${exp_unique_index[@]}" >> list_of_commands.txt
	done
done
