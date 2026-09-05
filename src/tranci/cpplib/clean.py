#!/usr/bin/env python3

import os
import glob

for name in ["src/main.x","src/main.exe"] + glob.glob("src/*.op"):
  if os.path.exists(name): os.remove(name) # remove this file
