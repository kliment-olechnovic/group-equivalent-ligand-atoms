#!/bin/bash

cd $(dirname "$0")

if [ ! -d "./testvenv" ]
then
	python3 -m venv testvenv

	source ./testvenv/bin/activate

	pip install rdkit

	pip install gemmi
fi

