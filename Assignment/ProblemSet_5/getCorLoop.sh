for FILE in ./SingleCluster_Data_L/L*.txt
do
    L=$(echo $FILE | sed 's/[^0-9]*//g')
    echo $L >> SingleCluster_Corr_L.txt
    ./get_tauA_M.out $FILE >> SingleCluster_Corr_L_tauA_M.txt
    ./get_tauA_E.out $FILE >> SingleCluster_Corr_L_tauA_E.txt
done
