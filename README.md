# LCS-papier
# for Linux
# You may need to build Krimp if that's the case go into Krimp/trunks and execute following commands line
./bootstrap.sh
make -Cbuild install

#if Res and experiments aren't directory please make them
# Then to execute the code from the LCS-papier directory
# c the number of transactions, n the number of sample that we will do and f the threshold 
python3 launch.py {database-name} {c} {n} {f}

# all results wanted would be in Res/{databsae-name}

# be carefull to save the directory with another name before relaunching another experiment
# the new experiment will overwrite every experiment on the same database
