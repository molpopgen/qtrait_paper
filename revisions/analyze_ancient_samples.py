import fwdpy11
import gzip
import sqlite3
import argparse
import sys
import os
import libsequence.variant_matrix as vm
import libsequence.summstats as sstats
import numpy as np
import pandas as pd

def make_parser():
    parser = argparse.ArgumentParser()
