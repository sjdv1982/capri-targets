python3 classify-entities.py on-hold-targets.json > classified-entities.json
python3 filter-targets-1.py classified-entities.json > filtered-targets-1.json 
python3 search-pdb.py filtered-targets-1.json search-pdb/ > searched-pdb.json
python3 analyze-searches.py searched-pdb.json > analyzed-searches.json
python3 classify-targets.py analyzed-searches.json > classified-targets.json
python3 filter-targets-2.py classified-targets.json > final-targets.json