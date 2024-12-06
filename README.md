# LCS-papier
# for Linux
# Building Krimp
You may need to build Krimp if that's the case go into Krimp/trunks and execute following commands line :
./bootstrap.sh
make -Cbuild install

# Running the code
if Res and experiments aren't directory please make them
Then to execute the code from the LCS-papier directory
c the number of transactions, n the number of sample that we will do and f the threshold  :

python3 launch.py {database-name} {c} {n} {f}
# Getting result
 all results wanted would be in Res/{database-name}

be carefull to save the directory with another name before relaunching another experiment
the new experiment will overwrite every experiment on the same database
