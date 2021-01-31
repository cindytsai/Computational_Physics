for FILE in ./HMC_lambda1E-6/*.txt
do
    #L=$(echo $FILE | sed 's/[^0-9]*//g')
    
    ./binning.out $FILE 1 >> HMC_Prop.txt
    ./binning.out $FILE 2 >> HMC_Sub.txt
done
