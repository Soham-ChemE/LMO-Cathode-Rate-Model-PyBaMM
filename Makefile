.PHONY: all install test clean

all:
	python scripts/run_all.py

install:
	pip install -r requirements.txt
	pip install -e .

test:
	python -m pytest -q

clean:
	rm -rf figures/*.png results/*.json data/processed/*.csv
