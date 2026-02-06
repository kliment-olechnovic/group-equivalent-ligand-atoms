#!/bin/bash

cd $(dirname "$0")

rm -rf "./output"
mkdir -p "./output"

####################################################################

./ensure_venv.bash

source ./testvenv/bin/activate

####################################################################

find ./input/ -type f \
| while read -r INFILE
do
	BASENAME="$(basename ${INFILE} | tr '.' '_')"
	
	python3 -B "../group-equivalent-ligand-atoms.py" "${INFILE}" > "./output/table_${BASENAME}.tsv"
done

################################################################################

git status -s ./output/

