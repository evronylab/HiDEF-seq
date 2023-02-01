#!/usr/bin/env python
import sys
from pbcore.io import IndexedBamReader
data = IndexedBamReader(sys.argv[1])
print(len(set(data.holeNumber)))
