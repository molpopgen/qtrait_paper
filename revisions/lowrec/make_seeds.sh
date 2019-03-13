#!/bin/bash

python3 ../makeseeds.py 256 lomu_seeds.txt $RANDOM
python3 ../makeseeds.py 256 midmu_seeds.txt $RANDOM
python3 ../makeseeds.py 256 himu_seeds.txt $RANDOM
