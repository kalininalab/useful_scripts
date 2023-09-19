list_file="$1"
dbpath="$2"
run_name="$3"

PanACota annotate -l "$list_file" -d "$dbpath" -r "./annotate_out/" -n "$run_name" --threads 32

PanACota pangenome -l "./annotate_out/LSTINFO-$list_file" -d "./annotate_out/Proteins/" -o "./pangenome_out/" -n "$run_name" --threads 32

PanACota corepers -p "./pangenome_out/PanGenome-$run_name.All.prt-clust-0.8-mode1.lst" -o "./corepers_out/"

PanACoTA align -c "./corepers_out/PersGenome_PanGenome-$run_name.All.prt-clust-0.8-mode1.lst-all_1.lst" -l "./annotate_out/LSTINFO-$list_file" -n "$run_name" -d "./annotate_out/" -o "./align_out" --threads 32

PanACoTA tree -a "./align_out/Phylo-$run_name/$run_name.nucl.grp.aln" -o "./tree/" --threads 32